#!/usr/bin/env python3
"""Conditional clean-IQ diagnostic, NOT a deployable/non-oracle estimator.

R_w(d) = mean_e rotation_e C_e[exp(i*omega*d) D_d primary_IQ].
Integer-shift responses commute with the linear frozen correlator, allowing
the exact finite windowed-sinc operator to be assembled without GPU reruns.
IQ16 re-quantization and CUDA floating-point reduction are not replicated;
agreement with independent GPU responses is therefore a mandatory gate.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
from dataclasses import dataclass, asdict
import json
import math
from pathlib import Path
import sys
import time

import numpy as np
from scipy.fft import fft, ifft, next_fast_len

import analyze_prn14_delay_mapping_sweep as search

HIGH = Path('<LOCAL_USER>/chh/GPU/High')
sys.path.insert(0, str(HIGH))
from audit_isolated_delay_profiles import load_responses, interp_uniform

MANIFEST = Path('<TMP_DIR>/prn14_bandlimited_interval_doulf9nl/manifest.json')
OUT = HIGH/'gps_l1ca_multisat_60s/isolated_delay_profile_audit/prn14_chain_response_control'


@dataclass(frozen=True)
class ResponseSpec:
    """Numerical definition of the receiver-response operator.

    A bank is valid only for the processing chain described by this object.
    The old script had these values as unrelated module constants, which made
    it too easy to reuse a bank with a different receiver configuration.
    """

    prn: int = 14
    sample_rate_hz: float = 25.0e6
    if_hz: float = -4.58e6
    chip_rate_hz: float = 1.023e6
    code_length_chips: int = 1023
    tap_first_chip: float = -2.0
    tap_step_chip: float = 0.01
    tap_count: int = 401
    delay_filter_half: int = 128
    delay_filter_beta: float = 10.0

    @property
    def chip_per_sample(self) -> float:
        return self.chip_rate_hz / self.sample_rate_hz

    @property
    def omega_per_sample(self) -> float:
        return 2.0 * np.pi * self.if_hz / self.sample_rate_hz

    @property
    def tap_grid(self) -> np.ndarray:
        return (self.tap_first_chip + self.tap_step_chip *
                np.arange(self.tap_count, dtype=np.float32))


DEFAULT_SPEC = ResponseSpec()
# Kept as compatibility aliases for the existing unit tests and analysis code.
H = DEFAULT_SPEC.chip_per_sample
OMEGA = DEFAULT_SPEC.omega_per_sample
WINDOWS = (0, 5, 8, 12, 13)
# Covers both columns over the unchanged tau/separation search domain,
# including all support of the 257-tap fractional-delay filter.
KMIN, KMAX = -158, 186


G2_DELAYS = (5, 6, 7, 8, 17, 18, 139, 140, 141, 251, 252, 254, 255,
             256, 257, 258, 469, 470, 471, 472, 473, 474, 509, 512, 513,
             514, 515, 516, 859, 860, 861, 862)


def ca_code(prn: int) -> np.ndarray:
    """GPS L1 C/A code in the same +/-1 orientation as the receiver source."""
    if not 1 <= prn <= len(G2_DELAYS):
        raise ValueError(f'GPS L1 C/A PRN must be in 1..32, got {prn}')
    a = [1]*10; b = [1]*10; g1 = []; g2 = []
    for _ in range(1023):
        g1.append(a[0]); g2.append(b[0])
        f1 = a[7]^a[0]; f2 = b[8]^b[7]^b[4]^b[2]^b[1]^b[0]
        a = a[1:]+[f1]; b = b[1:]+[f2]
    delay = 1023 - G2_DELAYS[prn - 1]
    bits = np.bitwise_xor(np.asarray(g1), np.roll(np.asarray(g2), -delay))
    return np.where(bits != 0, 1.0, -1.0)


def ca14():
    """Compatibility wrapper for the original PRN14-only diagnostic."""
    return ca_code(14)


def states_by_window(path, prn=14, windows=WINDOWS,
                     samples_per_window=100_000_000):
    groups = {}
    with Path(path).open() as f:
        for row in csv.DictReader(f):
            if int(row['prn']) != prn:
                continue
            w = int(row['input_start_sample'])//int(samples_per_window)
            if w not in windows:
                continue
            rows = groups.setdefault(w, [])
            if len(rows) < 100:
                rows.append(row)
    missing = [w for w in windows if len(groups.get(w, [])) != 100]
    if missing:
        raise ValueError(f'{path}: incomplete 100-epoch windows {missing}')
    return groups


def weights(row, taps, spec=DEFAULT_SPEC, fma=True):
    """CPU replica of float32 NCO arithmetic (not CUDA's approximate trig)."""
    n = np.arange(int(row['correlation_length_samples']), dtype=np.float32)
    step = np.float32(float(row['code_phase_step_chips']))
    rem = np.float32(float(row['rem_code_phase_chips']))
    # CUDA contracts multiplication+addition; subtraction is a second op.
    if fma:
        phase = (n[None,:].astype(float)*float(step)+taps[:,None].astype(float)).astype(np.float32)-rem
    else:
        phase = n[None,:]*step+taps[:,None]-rem
    phase = np.fmod(phase, np.float32(spec.code_length_chips))
    phase[phase < 0] += np.float32(spec.code_length_chips)
    # A tiny negative remainder can round to exactly 1023 after float32
    # addition. Use periodic wrapping here; the CUDA source has no guard.
    # Forward validation measures this and other CPU/GPU discrepancies.
    code = ca_code(spec.prn)[np.floor(phase).astype(int) % spec.code_length_chips]
    start = int(row['input_start_sample'])
    c0 = np.float32(math.fmod(float(row['rem_carrier_phase_rad'])+
                               spec.omega_per_sample * start, 2*np.pi))
    cs = np.float32(float(row['carrier_phase_step_rad'])+spec.omega_per_sample)
    cp = (n.astype(float)*float(cs)+float(c0)).astype(np.float32) if fma else n*cs+c0
    wipe = np.exp(-1j*cp.astype(float))
    return code*wipe[None,:]


def build_bank(window, out, fma=True, manifest_path=MANIFEST,
               spec=DEFAULT_SPEC, primary_path=None, state_path=None,
               samples_per_window=100_000_000):
    """Build one response bank from a calibration clean-IQ recording.

    The calibration recording is deliberately an explicit input.  It may not
    be the mixed test recording.  Reusing the test clean component is allowed
    only for the historical diagnostic and is labelled as such in metadata.
    """
    out = Path(out); manifest = json.loads(Path(manifest_path).read_text())
    raw_path = Path(primary_path or manifest['primary'])
    state_file = Path(state_path or manifest['state'])
    raw = np.memmap(raw_path, dtype='<i2', mode='r').reshape(-1,2)
    rows = states_by_window(
        state_file, prn=spec.prn, windows=(window,),
        samples_per_window=samples_per_window)[window]
    taps = spec.tap_grid
    bank = np.zeros((KMAX-KMIN+1, spec.tap_count), complex)
    started = time.monotonic()
    for epoch, row in enumerate(rows):
        start = int(row['input_start_sample']); count = int(row['correlation_length_samples'])
        if start - KMAX < 0 or start + count - KMIN > len(raw):
            raise ValueError(
                f'calibration IQ does not cover bank support for window {window}, '
                f'epoch {epoch}: start={start}, count={count}, samples={len(raw)}')
        part = raw[start-KMAX:start+count-KMIN].astype(float)
        signal = part[:,0]+1j*part[:,1]
        size = next_fast_len(len(signal)+count-1)
        spectrum = fft(signal, size)
        # A separate prompt determines the same primary-referenced rotation.
        wp = weights(row, taps[200:201], spec, fma)[0]
        p = np.dot(signal[KMAX:KMAX+count], wp)
        rotation = np.conj(p)/abs(p) if abs(p)>1e-6 else 1
        for first in range(0, spec.tap_count, 16):
            w = weights(row, taps[first:first+16], spec, fma)
            conv = ifft(fft(w[:,::-1],size,axis=1)*spectrum[None,:],axis=1)
            bank[:,first:first+16] += rotation*conv[:,count-1:count+KMAX-KMIN][:,::-1].T
        if (epoch+1)%25 == 0:
            print(f'bank window={window} epochs={epoch+1}/100 elapsed={time.monotonic()-started:.1f}s',flush=True)
    bank /= len(rows)
    out.mkdir(parents=True, exist_ok=True)
    metadata = {
        'model': 'receiver_chain_response_bank',
        'source_mode': 'clean_iq_calibration',
        'spec': asdict(spec),
        'kmin': KMIN, 'kmax': KMAX, 'fma': bool(fma),
        'manifest': str(manifest_path), 'calibration_primary': str(raw_path),
        'calibration_state': str(state_file), 'window': int(window),
        'samples_per_window': int(samples_per_window),
    }
    np.savez(out/f'bank_{window}.npz', bank=bank, kmin=KMIN, kmax=KMAX,
             fma=fma, metadata=json.dumps(metadata))
    return window


def _scalar_from_npz(data, name, default):
    if name not in data:
        return default
    value = np.asarray(data[name]).reshape(-1)[0]
    return value.item() if hasattr(value, 'item') else value


def load_response_model(path, expected_spec=None):
    """Load a bank and reject a calibration/test processing-chain mismatch."""
    with np.load(path, allow_pickle=False) as data:
        bank = np.asarray(data['bank'])
        kmin = int(_scalar_from_npz(data, 'kmin', KMIN))
        kmax = int(_scalar_from_npz(data, 'kmax', KMAX))
        metadata_raw = _scalar_from_npz(data, 'metadata', '{}')
    try:
        metadata = json.loads(str(metadata_raw))
    except json.JSONDecodeError as exc:
        raise ValueError(f'{path}: invalid bank metadata') from exc
    spec_data = metadata.get('spec', {})
    spec = ResponseSpec(**{
        key: spec_data[key] for key in asdict(DEFAULT_SPEC)
        if key in spec_data
    })
    if expected_spec is not None:
        for key in ('prn', 'sample_rate_hz', 'if_hz', 'chip_rate_hz',
                    'code_length_chips', 'tap_first_chip', 'tap_step_chip',
                    'tap_count', 'delay_filter_half', 'delay_filter_beta'):
            actual = getattr(spec, key)
            expected = getattr(expected_spec, key)
            if not np.isclose(actual, expected, rtol=0.0, atol=1.0e-12):
                raise ValueError(
                    f'{path}: bank {key}={actual} does not match expected {expected}')
    return ChainResponse(bank, kmin=kmin, kmax=kmax, spec=spec, metadata=metadata)


class ChainResponse:
    def __init__(self, bank, kmin=KMIN, kmax=KMAX, spec=DEFAULT_SPEC,
                 metadata=None):
        self.bank = bank
        self.kmin = int(kmin)
        self.kmax = int(kmax)
        self.spec = spec
        self.metadata = metadata or {}
        if bank is not None and bank.shape[0] != self.kmax - self.kmin + 1:
            raise ValueError('response bank support does not match kmin/kmax')
        self.table = None

    def exact(self, delay_chip):
        delay = np.atleast_1d(delay_chip)/self.spec.chip_per_sample
        integer = np.floor(delay).astype(int)
        frac = delay-integer
        half = self.spec.delay_filter_half
        r = np.arange(-half, half+1)
        h = (np.sinc(r[None,:]-frac[:,None]) *
             np.kaiser(2*half+1, self.spec.delay_filter_beta)[None,:])
        h /= h.sum(axis=1)[:,None]
        indices = integer[:,None]+r[None,:]-self.kmin
        if indices.min()<0 or indices.max()>=len(self.bank):
            raise ValueError('Candidate outside response-bank support')
        matrix = np.zeros((len(delay),len(self.bank)))
        np.put_along_axis(matrix,indices,h,axis=1)
        return (matrix@self.bank)*np.exp(1j*self.spec.omega_per_sample*delay)[:,None]

    def prepare(self):
        ticks = np.arange(-2400,4801)
        self.table = self.exact(ticks*.0005)

    def response(self, d):
        ticks = np.rint(np.asarray(d)/.0005).astype(int)
        if self.table is not None and np.allclose(np.asarray(d),ticks*.0005,atol=1e-10,rtol=0):
            if ticks.min() < -2400 or ticks.max() > 4800:
                raise ValueError('Candidate outside unchanged search domain')
            return self.table[ticks+2400]
        return self.exact(d)

    def profile(self, x, _primary, mixed, separation_grid, tau_grid):
        # Same four complex columns, LS implementation, and winner selection
        # as audit_isolated_delay_profiles.profile_for_window. ONLY b0/b1 differ.
        profile = np.full(len(separation_grid),np.inf)
        best_tau = np.full(len(separation_grid),np.nan)
        constant = np.broadcast_to(np.ones_like(x),(len(tau_grid),len(x)))
        linear = np.broadcast_to(x,(len(tau_grid),len(x)))
        b0 = self.response(-tau_grid)
        for idx,separation in enumerate(separation_grid):
            b1 = self.response(separation-tau_grid)
            design = np.stack((b0,b1,constant,linear),axis=2)
            # Scale each complex column before forming normal equations.  This
            # preserves the fitted subspace while preventing arbitrary signal
            # amplitude from dominating the numerical conditioning.
            scales = np.sqrt(np.einsum('nmc,nmc->nc', design.conj(), design).real)
            scales = np.maximum(scales, 1.0e-30)
            scaled = design / scales[:, None, :]
            gram = np.einsum('nmc,nmd->ncd',scaled.conj(),scaled)
            rhs = np.einsum('nmc,m->nc',scaled.conj(),mixed)
            try:
                scaled_coefficients = np.linalg.solve(gram,rhs[...,None])[...,0]
            except np.linalg.LinAlgError:
                scaled_coefficients = np.einsum(
                    'ncd,nd->nc', np.linalg.pinv(gram,rcond=1e-12), rhs)
            # A tiny residual can still be meaningless when the two response
            # columns are nearly indistinguishable.  Detect that situation in
            # the scaled design and recompute only those candidates with a
            # direct SVD least-squares solve.
            hermitian_gram = (gram + gram.conj().transpose(0, 2, 1)) * 0.5
            eigenvalues = np.linalg.eigvalsh(hermitian_gram)
            condition = np.sqrt(
                eigenvalues[:, -1] /
                np.maximum(eigenvalues[:, 0], 1.0e-14))
            unstable = (~np.isfinite(condition)) | (condition > 1.0e8)
            for candidate in np.flatnonzero(unstable):
                scaled_coefficients[candidate] = np.linalg.lstsq(
                    scaled[candidate], mixed, rcond=1.0e-12)[0]
            coefficients = scaled_coefficients / scales
            residual = mixed[None,:]-np.einsum('nmc,nc->nm',design,coefficients)
            costs = np.sum(np.abs(residual)**2,axis=1)
            winner = int(np.argmin(costs))
            profile[idx] = costs[winner]; best_tau[idx] = tau_grid[winner]
            self.last_diagnostics = {
                'max_scaled_condition': float(np.nanmax(condition)),
                'unstable_candidate_count': int(np.count_nonzero(unstable)),
                'tau_candidate_count': int(len(tau_grid)),
                'separation': float(separation),
            }
        return profile,best_tau,separation_grid.copy()


class SignalSimLinearChainResponse(ChainResponse):
    """Response bank using SignalSim's actual fractional-delay operator.

    ``ChainResponse`` is retained unchanged for reproducibility of the frozen
    sinc-based experiments.  SignalSim generates a delayed sample with an
    integer shift followed by linear interpolation, so new receiver-chain
    experiments must use this explicit model instead of silently reusing the
    old windowed-sinc operator.
    """

    def exact(self, delay_chip):
        delays = np.atleast_1d(np.asarray(delay_chip, dtype=float))
        delay_samples = delays / self.spec.chip_per_sample
        integer = np.floor(delay_samples).astype(int)
        fraction = delay_samples - integer
        index0 = integer - self.kmin
        index1 = integer + 1 - self.kmin
        if index0.min() < 0 or index1.max() >= len(self.bank):
            raise ValueError('Candidate outside response-bank support')
        response = ((1.0 - fraction)[:, None] * self.bank[index0] +
                    fraction[:, None] * self.bank[index1])
        return response * np.exp(
            1j * self.spec.omega_per_sample * delay_samples
        )[:, None]


def load_linear_response_model(path, expected_spec=None):
    """Load a bank with the SignalSim linear fractional-delay operator.

    The loader is intentionally separate from ``load_response_model`` so old
    frozen results continue to use their recorded operator explicitly.
    """
    model = load_response_model(path, expected_spec=expected_spec)
    return SignalSimLinearChainResponse(
        model.bank, kmin=model.kmin, kmax=model.kmax, spec=model.spec,
        metadata={**model.metadata, 'delay_operator': 'signalsim_linear'}
    )


def spec_from_manifest(path, prn=14):
    """Create a processing-chain spec from a simulator/receiver manifest."""
    data = json.loads(Path(path).read_text())
    sample_rate = float(data.get('sample_rate', DEFAULT_SPEC.sample_rate_hz))
    if_hz = float(data.get('if_hz', DEFAULT_SPEC.if_hz))
    chip_per_sample = data.get('chip_per_sample')
    chip_rate = (float(chip_per_sample) * sample_rate if chip_per_sample is not None
                 else DEFAULT_SPEC.chip_rate_hz)
    return ResponseSpec(prn=int(prn), sample_rate_hz=sample_rate,
                        if_hz=if_hz, chip_rate_hz=chip_rate)


def validate(out, manifest_path=MANIFEST, windows=WINDOWS, spec=DEFAULT_SPEC):
    manifest = json.loads(Path(manifest_path).read_text()); results = []
    models = {w:load_response_model(Path(out)/f'bank_{w}.npz', spec)
              for w in windows}
    for case in manifest['cases']:
        data = np.genfromtxt(case['response'],delimiter=',',names=True)
        for w in windows:
            a = data[data['window']==w]; x = a['delay_chip']
            p = a['primary_i']+1j*a['primary_q']; s = a['secondary_i']+1j*a['secondary_q']
            d = case['true_delay_chip']; pred = models[w].exact([d])[0]
            zero = models[w].exact([0])[0]
            shifted = interp_uniform(p,x[0],np.median(np.diff(x)),x+d)
            def gain_error(b):
                gain = np.vdot(b,s)/np.vdot(b,b)
                return float(np.linalg.norm(s-gain*b)/np.linalg.norm(s))
            results.append(dict(window=w,true_delay_chip=d,
                primary_relative_error=float(np.linalg.norm(p-zero)/np.linalg.norm(p)),
                secondary_relative_error=float(np.linalg.norm(s-pred)/np.linalg.norm(s)),
                chain_shape_error=gain_error(pred),shifted_template_shape_error=gain_error(shifted)))
    write_rows(out/'forward_validation.csv',results)
    summary = {k:max(r[k] for r in results) for k in ('primary_relative_error','secondary_relative_error','chain_shape_error','shifted_template_shape_error')}
    # Sub-percent gate is explicit; report actual errors, do not call bit exact.
    summary['passed'] = summary['primary_relative_error']<.001 and summary['secondary_relative_error']<.001
    (out/'forward_validation.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2),flush=True)
    if not summary['passed']:
        raise RuntimeError('CPU/GPU forward validation failed; do not interpret fitted delays')


def write_rows(path,rows):
    with path.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def fit_window(window, out, manifest_path=MANIFEST, spec=DEFAULT_SPEC):
    out=Path(out); model=load_response_model(out/f'bank_{window}.npz', spec);model.prepare()
    # Reuse the existing coarse-to-fine driver verbatim, no truth initialization.
    search.profile_for_window=model.profile
    result=[]
    for case in json.loads(Path(manifest_path).read_text())['cases']:
        path=out/f'fit_{window}_{case["true_delay_chip"]:.9f}.json'
        if path.exists():
            row=json.loads(path.read_text())
        else:
            row=search.evaluate_case(case['true_delay_chip'],window,case['response'])
            path.write_text(json.dumps(row,indent=2)+'\n')
        result.append(row)
        print(json.dumps(row),flush=True)
    return result


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('stage',choices=['bank','validate','fit'])
    parser.add_argument('--output-dir',type=Path,default=OUT)
    parser.add_argument('--manifest',type=Path,default=MANIFEST,
                        help='independent mixed-test manifest')
    parser.add_argument('--calibration-manifest',type=Path,
                        help='clean calibration manifest used only for bank stage')
    parser.add_argument('--calibration-primary',type=Path,
                        help='override clean calibration IQ path')
    parser.add_argument('--calibration-state',type=Path,
                        help='override clean calibration state CSV path')
    parser.add_argument('--allow-same-calibration',action='store_true',
                        help='explicitly allow the historical same-recording diagnostic')
    parser.add_argument('--prn',type=int,default=14)
    parser.add_argument('--workers',type=int,default=5)
    parser.add_argument('--windows',type=int,nargs='+',default=list(WINDOWS))
    parser.add_argument('--no-fma',action='store_true')
    args=parser.parse_args();out=args.output_dir;out.mkdir(parents=True,exist_ok=True)
    spec = spec_from_manifest(args.manifest, args.prn)
    windows = tuple(args.windows)
    if args.stage=='validate':
        validate(out, args.manifest, windows, spec);return
    if args.stage=='fit':
        validate(out, args.manifest, windows, spec)
    calibration_manifest = args.calibration_manifest or args.manifest
    calibration_spec = spec_from_manifest(calibration_manifest, args.prn)
    for field in ('sample_rate_hz', 'if_hz', 'chip_rate_hz'):
        if not np.isclose(getattr(spec, field), getattr(calibration_spec, field),
                          rtol=0.0, atol=1.0e-9):
            raise SystemExit(
                f'calibration/test processing-chain mismatch in {field}: '
                f'{getattr(calibration_spec, field)} vs {getattr(spec, field)}')
    if args.stage == 'bank' and not args.allow_same_calibration:
        test_data = json.loads(args.manifest.read_text())
        calibration_data = json.loads(calibration_manifest.read_text())
        test_primary = Path(test_data['primary']).resolve()
        calibration_primary = Path(
            args.calibration_primary or calibration_data['primary']).resolve()
        if calibration_primary == test_primary:
            raise SystemExit(
                'refusing to build a deployable bank from the test clean-IQ file; '
                'provide an independent --calibration-manifest/--calibration-primary '
                'or explicitly pass --allow-same-calibration for the historical diagnostic')
    result=[]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        jobs=[pool.submit(
            build_bank, w, out, not args.no_fma, calibration_manifest, spec,
            args.calibration_primary, args.calibration_state)
            if args.stage=='bank' else pool.submit(
                fit_window, w, out, args.manifest, spec)
            for w in windows]
        for job in as_completed(jobs):
            r=job.result()
            if args.stage=='fit':result.extend(r)
    if args.stage=='fit':
        result.sort(key=lambda r:(r['true_delay_chip'],r['window']))
        write_rows(out/'chain_global_profile.csv',result)
        errors=np.array([r['global_bias_chip'] for r in result])
        summary=dict(scope='independent-calibration response model',count=len(result),rmse_chip=float(np.sqrt(np.mean(errors**2))),max_abs_error_chip=float(max(abs(errors))),within_0p01=int(np.sum(abs(errors)<=.01000000001)),spec=asdict(spec))
        (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
        print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
