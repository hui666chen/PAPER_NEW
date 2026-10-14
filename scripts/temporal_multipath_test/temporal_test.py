#!/usr/bin/env python3
"""Short-delay multipath: Hatch headroom + geometry-model estimator (correlator level, no IQ).

Physics.  A reflector at height h below the antenna (horizontal surface) gives
    delay(t) = 2 h sin(el(t))            [m]
    phase(t) = 2 pi delay(t) / lambda + phi0
so the relative phase rotates as the satellite moves.  A single 1-s curve
cannot separate a short-delay reflection from a shift of the direct path, but
a few minutes of curves can, because (h, alpha, phi0) are constants while
el(t) is known from the ephemeris.

Per scenario (satellite pass x h x alpha) this script simulates 1-s 401-tap
curves with the calibrated receiver response, then compares
    EML             0.1-chip power early-minus-late, no smoothing
    HATCH           EML + Hatch carrier smoothing (window --hatch, default 100 s)
    MODEL           EML minus the bias predicted by a causal fit of (h, alpha, phi0)
                    on the previous --window seconds of curves
    MODEL_HATCH     MODEL code + model-corrected carrier, then the same Hatch filter

Frozen decision (protocols/TEMPORAL_SHORT_DELAY_TEST_PROTOCOL.md) is evaluated
at the end and written to TEMPORAL_REPORT.md.

Run on the experiment PC (real response bank, real ephemeris):
    python3 scripts/temporal_multipath_test/temporal_test.py --workers 4
Self-test (synthetic response and orbits, ~3 min):
    python3 scripts/temporal_multipath_test/temporal_test.py --self-test
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
import json
from pathlib import Path
import sys
import time

import numpy as np
from scipy.optimize import minimize

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(HERE))
from orbit import synthetic_trajectories, trajectories  # noqa: E402

C_LIGHT = 299_792_458.0
CHIP_M = C_LIGHT / 1.023e6
LAMBDA = C_LIGHT / 1575.42e6
DEFAULT_PAPER = Path("/home/zoxoai/Documents/ChatGPT/paper")
DEFAULT_NAV = Path("/home/zoxoai/chh/GNSS/SignalSim-main/EphData/BRDC00IGS_R_20211700000_01D_MN.rnx")
SCENE = {"start": "2021-06-19T10:05:30", "lat": 37.352721, "lon": -121.915773, "alt": 20.0}
EML_SPACING = 0.10
S_GRID = np.round(np.arange(-0.30, 0.8000001, 0.001), 6)
D_GRID = np.round(np.arange(0.0, 0.4000001, 0.002), 6)
A_GRID = np.round(np.arange(0.0, 0.7000001, 0.05), 6)
P_GRID = np.arange(64) * (2 * np.pi / 64)
SEED = 2026101430
PRIMARY_H = (3.0, 5.0, 8.0, 12.0)
PRIMARY_A = (0.3, 0.5)


# ------------------------------------------------------------------ response models
class TriangleModel:
    """Band-limited triangle (self-test only)."""

    def __init__(self, tap):
        self.tap = tap
        self.k = np.linspace(-0.04, 0.04, 9)

    def exact(self, d):
        d = np.atleast_1d(np.asarray(d, float))
        x = self.tap[None, :, None] + d[:, None, None] + self.k[None, None, :]
        return np.maximum(0.0, 1.0 - np.abs(x)).mean(axis=2).astype(np.complex128)


def load_model(args):
    if args.self_test:
        tap = -2.0 + 0.01 * np.arange(401)
        return TriangleModel(tap), tap
    for p in (args.paper, args.paper / "PAPER_NEW" / "scripts", args.paper / "gps_l1_prn09_shared_ratio_v2",
              REPO / "scripts"):
        sys.path.insert(0, str(p))
    try:
        from receiver_chain_response_control import ResponseSpec, load_response_model
    except ImportError as exc:
        raise SystemExit(f"cannot import receiver_chain_response_control ({exc}); run on the experiment PC")
    spec = ResponseSpec(prn=9)
    return load_response_model(args.bank, expected_spec=spec), spec.tap_grid.astype(np.float64)


# ------------------------------------------------------------------ vectorised EML
def eml_roots(curves: np.ndarray, tap: np.ndarray, d: float = EML_SPACING) -> np.ndarray:
    """Power EML root nearest zero, same definition as correlator_baselines.eml_estimate."""
    grid = np.arange(-0.5, 0.5 + 0.001, 0.002)
    p = np.abs(np.atleast_2d(curves)) ** 2
    t0, dt = float(tap[0]), float(tap[1] - tap[0])

    def interp(x):
        f = (x - t0) / dt
        i = np.clip(np.floor(f).astype(int), 0, len(tap) - 2)
        w = f - i
        return p[:, i] * (1 - w) + p[:, i + 1] * w

    v = interp(grid - d / 2) - interp(grid + d / 2)
    a, b = v[:, :-1], v[:, 1:]
    cross = (a * b < 0) | (a == 0)
    loc = grid[:-1] + np.where(a == 0, 0.0, np.abs(a) / (np.abs(a) + np.abs(b) + 1e-300)) * (grid[1] - grid[0])
    dist = np.where(cross, np.abs(loc), np.inf)
    k = np.argmin(dist, axis=1)
    out = loc[np.arange(len(p)), k]
    none = ~np.isfinite(dist[np.arange(len(p)), k])
    if none.any():
        out[none] = grid[np.argmin(np.abs(v[none]), axis=1)]
    return out


# ------------------------------------------------------------------ shared tables
class Tables:
    def __init__(self, model, tap):
        self.tap = tap
        one = np.ones(len(tap), complex)
        self.q, _ = np.linalg.qr(np.column_stack((one, tap.astype(complex))), mode="reduced")
        r0 = np.asarray(model.exact([0.0]))[0]
        self.k0 = int(np.argmin(np.abs(tap)))
        self.scale = 1.0 / np.abs(r0).max()
        self.r0 = r0 * self.scale                                  # direct, unit peak
        self.rd = np.asarray(model.exact(D_GRID)) * self.scale     # reflections on D_GRID
        self.c0 = float(eml_roots(self.r0[None], tap)[0])
        ts = np.asarray(model.exact(S_GRID)).T * self.scale
        self.ts = ts - self.q @ (self.q.conj().T @ ts)             # projected templates (taps x S)
        p0 = self.r0 - self.q @ (self.q.conj().T @ self.r0)
        pd = (self.rd.T - self.q @ (self.q.conj().T @ self.rd.T))
        self.n0 = float(np.vdot(p0, p0).real)
        self.nd = np.einsum("ij,ij->j", pd.conj(), pd).real
        self.x0 = p0.conj() @ pd                                    # <T(0), T(delta)>
        # EML bias table F[delta, alpha, phase] (chips) on noiseless composites
        f = np.empty((len(D_GRID), len(A_GRID), len(P_GRID)))
        for ia, a in enumerate(A_GRID):
            rho = a * np.exp(1j * P_GRID)
            comp = self.r0[None, None, :] + rho[None, :, None] * self.rd[:, None, :]
            f[:, ia, :] = (eml_roots(comp.reshape(-1, len(tap)), tap) - self.c0).reshape(len(D_GRID), len(P_GRID))
        self.f = f
        self.prompt_d = self.rd[:, self.k0] / self.r0[self.k0]     # reflection/direct at prompt

    def bias(self, delta, alpha, phase):
        """Trilinear interpolation of F (chips); delta in chips, phase in rad."""
        x = np.clip(np.asarray(delta) / 0.002, 0, len(D_GRID) - 1.0001)
        y = np.clip(np.asarray(alpha) / 0.05, 0, len(A_GRID) - 1.0001)
        z = (np.asarray(phase) % (2 * np.pi)) / (2 * np.pi) * len(P_GRID)
        i, j, k = x.astype(int), y.astype(int), z.astype(int) % len(P_GRID)
        fx, fy, fz = x - i, y - j, z - np.floor(z)
        k1 = (k + 1) % len(P_GRID)
        f = self.f
        c00 = f[i, j, k] * (1 - fx) + f[i + 1, j, k] * fx
        c01 = f[i, j, k1] * (1 - fx) + f[i + 1, j, k1] * fx
        c10 = f[i, j + 1, k] * (1 - fx) + f[i + 1, j + 1, k] * fx
        c11 = f[i, j + 1, k1] * (1 - fx) + f[i + 1, j + 1, k1] * fx
        return (c00 * (1 - fz) + c01 * fz) * (1 - fy) + (c10 * (1 - fz) + c11 * fz) * fy

    def lin(self, table, delta):
        x = np.clip(np.asarray(delta) / 0.002, 0, len(D_GRID) - 1.0001)
        i = x.astype(int)
        w = x - i
        return table[i] * (1 - w) + table[i + 1] * w

    def carrier_err_m(self, delta, alpha, phase):
        ratio = 1.0 + alpha * np.exp(1j * phase) * self.lin(self.prompt_d, delta)
        return np.angle(ratio) / (2 * np.pi) * LAMBDA


def tap_noise_factor(tap):
    """Factor of correlator noise across taps: covariance = chip autocorrelation (triangle)."""
    dx = np.abs(tap[:, None] - tap[None, :])
    cov = np.maximum(0.0, 1.0 - dx)
    w, v = np.linalg.eigh(cov)
    return v * np.sqrt(np.clip(w, 0.0, None))[None, :]


# ------------------------------------------------------------------ estimator
class GeometryFit:
    """Causal fit of (h, alpha, phi0) on a window of curves."""

    H = np.arange(1.0, 30.0001, 0.5)
    A = np.arange(0.05, 0.6001, 0.05)
    P = np.arange(24) * (2 * np.pi / 24)

    def __init__(self, tab: Tables, corr: np.ndarray, e_hat: np.ndarray, sin_el: np.ndarray):
        self.tab, self.corr, self.e_hat, self.sin_el = tab, corr, e_hat, sin_el

    def _cinterp(self, t_idx, s):
        f = (s - S_GRID[0]) / 0.001
        i = np.clip(np.floor(f).astype(int), 0, len(S_GRID) - 2)
        w = f - i
        return self.corr[t_idx, i] * (1 - w) + self.corr[t_idx, i + 1] * w

    def cost(self, h, a, p0, idx):
        """h, a, p0 broadcastable to (B, 1); idx (W,) epochs.  Returns (B,)."""
        tab = self.tab
        dm = 2.0 * h * self.sin_el[idx][None, :]
        delta = dm / CHIP_M
        phase = 2 * np.pi * dm / LAMBDA + p0
        rho = a * np.exp(1j * phase)
        tau = self.e_hat[idx][None, :] - tab.c0 - tab.bias(delta, np.broadcast_to(a, delta.shape), phase)
        tt = np.broadcast_to(idx[None, :], delta.shape)
        num = self._cinterp(tt, -tau) + np.conj(rho) * self._cinterp(tt, delta - tau)
        den = tab.n0 + 2 * np.real(rho * tab.lin(tab.x0, delta)) + (np.abs(rho) ** 2) * tab.lin(tab.nd, delta)
        return -(np.abs(num).sum(axis=1)) ** 2 / den.sum(axis=1)

    def cost_null(self, idx):
        tau = self.e_hat[idx] - self.tab.c0
        num = self._cinterp(idx, -tau)
        return float(-(np.abs(num).sum()) ** 2 / (len(idx) * self.tab.n0))

    def fit(self, idx):
        hh, aa, pp = np.meshgrid(self.H, self.A, self.P, indexing="ij")
        hh, aa, pp = hh.ravel()[:, None], aa.ravel()[:, None], pp.ravel()[:, None]
        best, best_j = None, np.inf
        for s in range(0, len(hh), 2048):
            j = self.cost(hh[s:s + 2048], aa[s:s + 2048], pp[s:s + 2048], idx)
            k = int(np.argmin(j))
            if j[k] < best_j:
                best_j, best = float(j[k]), (float(hh[s + k, 0]), float(aa[s + k, 0]), float(pp[s + k, 0]))

        def f(x):
            h, a, p = x
            if not (0.5 <= h <= 40 and 0.0 <= a <= 0.69):
                return 1e300
            return float(self.cost(np.array([[h]]), np.array([[a]]), np.array([[p]]), idx)[0])

        res = minimize(f, np.array(best), method="Nelder-Mead",
                       options={"xatol": 1e-3, "fatol": 1e-9 * abs(best_j), "maxiter": 400})
        if res.fun < best_j:
            best, best_j = tuple(res.x), float(res.fun)
        j0 = self.cost_null(idx)
        return best, best_j, j0


# ------------------------------------------------------------------ simulation
def hatch(code, carrier, n):
    s = np.empty_like(code)
    s[0] = code[0]
    for k in range(1, len(code)):
        w = min(k + 1, n)
        s[k] = code[k] / w + (w - 1) / w * (s[k - 1] + carrier[k] - carrier[k - 1])
    return s


def simulate(job):
    t0 = time.time()
    tab, chol, sigma = job["tab"], job["chol"], job["sigma"]
    rng = np.random.default_rng(job["seed"])
    el = job["el"]
    sin_el = np.sin(np.deg2rad(el))
    n = len(el)
    h, a, p0 = job["h"], job["alpha"], job["phi0"]
    dm = 2 * h * sin_el
    delta, phase = dm / CHIP_M, 2 * np.pi * dm / LAMBDA + p0
    # noiseless composite curves
    xi = np.clip(delta / 0.002, 0, len(D_GRID) - 1.0001)
    i = xi.astype(int)
    w = (xi - i)[:, None]
    rd_t = tab.rd[i] * (1 - w) + tab.rd[i + 1] * w
    curves = tab.r0[None, :] + (a * np.exp(1j * phase))[:, None] * rd_t
    noise = (rng.normal(size=(n, len(tab.tap))) + 1j * rng.normal(size=(n, len(tab.tap)))) / np.sqrt(2)
    curves = curves + sigma * noise @ chol.T
    if job["mismatch"] is not None:
        shapes = job["mismatch"]
        pick = rng.integers(len(shapes), size=n // 60 + 1)
        curves = curves + shapes[np.repeat(pick, 60)[:n]]
    # common gain/phase rotation (receiver phase), irrelevant to every method
    curves *= np.exp(1j * rng.uniform(0, 2 * np.pi))
    e_hat = eml_roots(curves, tab.tap)
    true_code_err = (tab.bias(delta, np.full(n, a), phase)) * CHIP_M
    carrier_err = tab.carrier_err_m(delta, a, phase) + rng.normal(0, 0.002, n)
    eml_err = (e_hat - tab.c0) * CHIP_M
    hat_err = hatch(eml_err, carrier_err, job["hatch"])
    # model estimator (causal)
    yperp = curves.T - tab.q @ (tab.q.conj().T @ curves.T)
    corr = (tab.ts.conj().T @ yperp).T          # (n, S)
    fitter = GeometryFit(tab, corr, e_hat, sin_el)
    win, every = job["window"], job["refit"]
    model_bias = np.zeros(n)
    model_carrier = np.zeros(n)
    fits = []
    for k in range(win, n, every):
        idx = np.arange(k - win, k)
        (hh, aa, pp), j1, j0 = fitter.fit(idx)
        gain = (j0 - j1) / abs(j0)
        use = gain > job["gate"]
        seg = slice(k, min(k + every, n))
        if use:
            dms = 2 * hh * sin_el[seg]
            ds, ps = dms / CHIP_M, 2 * np.pi * dms / LAMBDA + pp
            model_bias[seg] = tab.bias(ds, np.full(len(ds), aa), ps) * CHIP_M
            model_carrier[seg] = tab.carrier_err_m(ds, aa, ps)
        fits.append({"t": k, "h": hh, "alpha": aa, "phi0": pp, "rel_gain": gain, "used": bool(use)})
    mod_err = eml_err - model_bias
    modh_err = hatch(mod_err, carrier_err - model_carrier, job["hatch"])
    m = np.arange(n) >= win + job["hatch"]
    rms = lambda x: float(np.sqrt(np.mean(x[m] ** 2)))
    used = [f for f in fits if f["used"]]
    out = {k: job[k] for k in ("prn", "h", "alpha", "phi0")}
    out.update({
        "delay_chip_min": float(delta.min()), "delay_chip_max": float(delta.max()),
        "fade_period_s": float(np.median(2 * np.pi / np.maximum(np.abs(np.gradient(phase)), 1e-12))),
        "rms_true_mp_m": rms(true_code_err), "rms_eml_m": rms(eml_err), "rms_hatch_m": rms(hat_err),
        "rms_model_m": rms(mod_err), "rms_model_hatch_m": rms(modh_err),
        "fit_used_fraction": len(used) / max(len(fits), 1),
        "h_abs_err_median_m": float(np.median([abs(f["h"] - h) for f in used])) if used and a > 0 else None,
        "alpha_abs_err_median": float(np.median([abs(f["alpha"] - a) for f in used])) if used and a > 0 else None,
        "rel_gains": [f["rel_gain"] for f in fits],
        "seconds": round(time.time() - t0, 1),
    })
    return out


def calibrate_sigma(tab, chol, target_m, rng, n=400):
    noise = (rng.normal(size=(n, len(tab.tap))) + 1j * rng.normal(size=(n, len(tab.tap)))) / np.sqrt(2)
    e = eml_roots(tab.r0[None, :] + noise @ chol.T * 1e-3, tab.tap)
    return target_m / (np.std(e) * CHIP_M) * 1e-3


def mismatch_shapes(args, tab):
    """Relative residual shapes of real clean 100-ms windows after a single-path fit."""
    if args.self_test or args.no_mismatch:
        return None
    sys.path.insert(0, str(REPO / "scripts" / "pvt_information_test"))
    from common import read_epoch_matrix
    from stage_a_cases import BANKS, CLEAN_WINDOWS
    cal = args.calib_dir
    shapes = []
    for cell in BANKS:
        path = cal / "open_loop" / cell / "prn09_epoch_responses.csv"
        if not path.exists():
            print(f"[warn] {path} missing; skipping mismatch shapes for {cell}")
            continue
        for wnd, ep in read_epoch_matrix(path, CLEAN_WINDOWS[cell]).items():
            y = ep.mean(axis=0)
            yp = y - tab.q @ (tab.q.conj().T @ y)
            best = None
            for s in np.arange(-0.1, 0.1001, 0.001):
                k = int(round((s - S_GRID[0]) / 0.001))
                t = tab.ts[:, k]
                g = np.vdot(t, yp) / np.vdot(t, t)
                r = yp - g * t
                c = float(np.vdot(r, r).real)
                if best is None or c < best[0]:
                    best = (c, r / g)
            shapes.append(best[1])
    if not shapes:
        return None
    print(f"[info] {len(shapes)} real mismatch/noise shapes loaded (100-ms windows, unscaled = conservative)")
    return np.asarray(shapes)


def decide(rows):
    prim = [r for r in rows if r["h"] in PRIMARY_H and r["alpha"] in PRIMARY_A]
    clean = [r for r in rows if r["alpha"] == 0]
    agg = lambda rs, k: float(np.sqrt(np.mean([r[k] ** 2 for r in rs]))) if rs else float("nan")
    hatch_p, modh_p = agg(prim, "rms_hatch_m"), agg(prim, "rms_model_hatch_m")
    indiv = np.mean([r["rms_model_hatch_m"] <= 0.8 * r["rms_hatch_m"] for r in prim]) if prim else 0.0
    clean_h, clean_m = agg(clean, "rms_hatch_m"), agg(clean, "rms_model_hatch_m")
    clean_ratio = clean_m / clean_h if clean else float("nan")
    clean_ok = (clean_m <= max(1.05 * clean_h, clean_h + 0.05)) if clean else False
    if hatch_p < 0.8:
        dec = "STOP_NO_HEADROOM"
    elif modh_p <= 0.6 * hatch_p and indiv >= 0.75 and clean_ok:
        dec = "PASS_GO_TO_IQ"
    else:
        dec = "STOP_METHOD_NOT_GOOD_ENOUGH"
    return {"decision": dec, "primary_hatch_rms_m": hatch_p, "primary_model_hatch_rms_m": modh_p,
            "primary_ratio": modh_p / hatch_p if hatch_p else float("nan"),
            "primary_individual_pass_fraction": float(indiv), "clean_model_hatch_over_hatch": clean_ratio,
            "rule": "Hatch RMS (primary) >= 0.8 m; MODEL_HATCH <= 0.6 x HATCH aggregate, >= 75% of primary "
                    "scenarios <= 0.8 x HATCH, and clean MODEL_HATCH <= HATCH + max(5%, 0.05 m)",
            "clean_hatch_rms_m": clean_h, "clean_model_hatch_rms_m": clean_m}


def signalsim_segments(prn, el, start_s, h, a, phi0_deg, step=1):
    seg = []
    for k in range(0, len(el), step):
        dm = 2 * h * np.sin(np.deg2rad(el[k]))
        seg.append({"system": "GPS", "signal": "L1CA", "svid": int(prn), "start": float(start_s + k),
                    "end": float(start_s + k + step), "delayChip": round(dm / CHIP_M, 6), "amplitude": a,
                    "phaseDeg": round((360.0 * dm / LAMBDA + phi0_deg) % 360.0, 4)})
    return seg


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--paper", type=Path, default=DEFAULT_PAPER)
    ap.add_argument("--bank", type=Path, default=None)
    ap.add_argument("--calib-dir", type=Path, default=None)
    ap.add_argument("--nav", type=Path, default=DEFAULT_NAV)
    ap.add_argument("--start", default=SCENE["start"])
    ap.add_argument("--duration", type=int, default=1800)
    ap.add_argument("--min-el", type=float, default=15.0)
    ap.add_argument("--max-sats", type=int, default=4)
    ap.add_argument("--heights", type=float, nargs="*", default=[3.0, 5.0, 8.0, 12.0, 20.0])
    ap.add_argument("--alphas", type=float, nargs="*", default=[0.0, 0.3, 0.5])
    ap.add_argument("--phi0-deg", type=float, default=180.0)
    ap.add_argument("--hatch", type=int, default=100)
    ap.add_argument("--window", type=int, default=120)
    ap.add_argument("--refit", type=int, default=30)
    ap.add_argument("--eml-noise-1s", type=float, default=0.5)
    ap.add_argument("--no-mismatch", action="store_true")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--out", type=Path, default=REPO / "results" / "temporal_multipath_test")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    cal = args.paper / "gps_l1_prn09_short_delay_pilot_20261013" / "calibration_v2_after_acquisition_amendment"
    args.bank = args.bank or cal / "response_banks" / "bank_39.npz"
    args.calib_dir = args.calib_dir or cal
    out = args.out / ("self_test" if args.self_test else "")
    out.mkdir(parents=True, exist_ok=True)
    if args.self_test:
        args.duration, args.max_sats, args.heights, args.alphas = 900, 2, [5.0, 12.0], [0.0, 0.5]
    rng = np.random.default_rng(SEED)
    model, tap = load_model(args)
    print("[info] building tables ...", flush=True)
    tab = Tables(model, tap)
    chol = tap_noise_factor(tap)
    sigma = calibrate_sigma(tab, chol, args.eml_noise_1s, rng)
    mism = mismatch_shapes(args, tab)
    if args.self_test:
        sats = synthetic_trajectories(args.duration, args.min_el, args.max_sats)
        geometry = "synthetic passes"
    else:
        start = datetime.fromisoformat(args.start).replace(tzinfo=timezone.utc)
        sats = trajectories(args.nav, start, args.duration, SCENE["lat"], SCENE["lon"], SCENE["alt"], args.min_el)
        geometry = f"{args.nav.name} from {args.start} UTC"
    prns = sorted(sats, key=lambda p: -float(np.mean(sats[p][1])))[: args.max_sats]
    if not prns:
        raise SystemExit("no satellite above --min-el for the whole span")
    # Clean-only gate calibration: 99th percentile of the relative cost gain on alpha = 0.
    base = {"tab": tab, "chol": chol, "sigma": sigma, "mismatch": mism, "phi0": np.deg2rad(args.phi0_deg),
            "hatch": args.hatch, "window": args.window, "refit": args.refit}
    jobs_clean = [dict(base, prn=p, el=sats[p][1], h=5.0, alpha=0.0, gate=-np.inf, seed=SEED + 1000 + i)
                  for i, p in enumerate(prns)]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        cal_rows = list(pool.map(simulate, jobs_clean))
    gains = np.concatenate([r["rel_gains"] for r in cal_rows])
    gate = float(np.quantile(gains, 0.99, method="higher"))
    print(f"[info] clean-calibrated gate (99th pct relative cost gain) = {gate:.4g}", flush=True)
    jobs = []
    for i, p in enumerate(prns):
        for h in args.heights:
            for a in args.alphas:
                jobs.append(dict(base, prn=p, el=sats[p][1], h=h, alpha=a, gate=gate,
                                 seed=SEED + 10_000 + len(jobs)))
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        rows = list(pool.map(simulate, jobs))
    for r in rows:
        print(f"PRN{r['prn']:>3} h={r['h']:>4} a={r['alpha']:.1f} delay {r['delay_chip_min']:.3f}-"
              f"{r['delay_chip_max']:.3f} chip fade {r['fade_period_s']:.0f}s | EML {r['rms_eml_m']:.2f} "
              f"HATCH {r['rms_hatch_m']:.2f} MODEL {r['rms_model_m']:.2f} MODEL_HATCH {r['rms_model_hatch_m']:.2f}")
    dec = decide(rows)
    for r in rows:
        r.pop("rel_gains", None)
    payload = {"status": "complete", "self_test": args.self_test, "iq_generated": False, "geometry": geometry,
               "satellites": prns, "noise": {"eml_1s_m": args.eml_noise_1s, "mismatch_shapes": 0 if mism is None else len(mism)},
               "gate": gate, "settings": {k: (str(v) if isinstance(v, Path) else v) for k, v in vars(args).items()},
               "decision": dec, "rows": rows}
    (out / "temporal_summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    lines = ["# Short-delay temporal test (correlator level, no IQ)", "", f"Decision: **{dec['decision']}**  ",
             f"Rule: {dec['rule']}", "",
             f"- Primary (h {PRIMARY_H}, alpha {PRIMARY_A}): HATCH {dec['primary_hatch_rms_m']:.3f} m, "
             f"MODEL_HATCH {dec['primary_model_hatch_rms_m']:.3f} m, ratio {dec['primary_ratio']:.3f}, "
             f"individual pass {dec['primary_individual_pass_fraction']:.2f}",
             f"- Clean MODEL_HATCH / HATCH: {dec['clean_model_hatch_over_hatch']:.3f}",
             f"- Geometry: {geometry}; satellites {prns}; EML 1-s noise {args.eml_noise_1s} m; "
             f"real mismatch shapes: {payload['noise']['mismatch_shapes']}; gate {gate:.4g}", "",
             "| PRN | h (m) | alpha | delay (chip) | fade (s) | true MP | EML | HATCH | MODEL | MODEL_HATCH | fits used | median |h err| |",
             "|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for r in rows:
        he = "" if r["h_abs_err_median_m"] is None else f"{r['h_abs_err_median_m']:.2f}"
        lines.append(f"| {r['prn']} | {r['h']} | {r['alpha']} | {r['delay_chip_min']:.3f}-{r['delay_chip_max']:.3f} | "
                     f"{r['fade_period_s']:.0f} | {r['rms_true_mp_m']:.2f} | {r['rms_eml_m']:.2f} | {r['rms_hatch_m']:.2f} | "
                     f"{r['rms_model_m']:.2f} | {r['rms_model_hatch_m']:.2f} | {r['fit_used_fraction']:.2f} | {he} |")
    lines += ["", "RMS values in metres, after the first window + Hatch warm-up. Curves are simulated with the "
              "calibrated response (the same one the estimator uses), so this is an upper bound on what the "
              "method can do; real-residual shapes are added when available to include response mismatch."]
    (out / "TEMPORAL_REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    if dec["decision"] == "PASS_GO_TO_IQ" and not args.self_test:
        best = min((r for r in rows if r["alpha"] > 0), key=lambda r: r["rms_model_hatch_m"] / r["rms_hatch_m"])
        seg = signalsim_segments(best["prn"], sats[best["prn"]][1], 0.0, best["h"], best["alpha"], args.phi0_deg)
        (out / f"signalsim_multipath_PRN{best['prn']}_h{best['h']}_a{best['alpha']}.json").write_text(
            json.dumps({"multipath": seg}, indent=1), encoding="utf-8")
    print(f"[decision] {dec['decision']} -> {out}")


if __name__ == "__main__":
    main()
