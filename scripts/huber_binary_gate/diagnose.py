"""Post-hoc diagnostic only; does not alter fits, gates, or PVT outputs."""
import json
import subprocess
from pathlib import Path
import numpy as np
from run import ROOT, BASE, SCENES
from score import helper, OUT

DEST = ROOT/'results/huber_correction_diagnostic_20261021'


def observations(path):
    records = {}
    epoch = None
    for line in path.read_text().splitlines():
        if line.startswith('>'):
            epoch = helper.parse_rinex_epoch(line)
            second = helper.signal_relative_second_utc(epoch-helper.timedelta(seconds=18))
            records[second] = dict(epoch=epoch, sats={})
        elif epoch is not None and line.startswith('G'):
            fields = [line[3+16*i:17+16*i].strip() for i in range(4)]
            if fields[0]:
                records[second]['sats'][line[:3]] = [float(x) if x else None for x in fields]
    return records


def summary(x):
    x = np.asarray(x, float)
    return dict(n=len(x), mean=float(np.mean(x)), median=float(np.median(x)),
                std=float(np.std(x)), min=float(np.min(x)), max=float(np.max(x)))


def main():
    DEST.mkdir(exist_ok=True)
    seal = json.loads((OUT/'PREDICTIONS_SEAL.json').read_text())
    clean_path = Path(seal['predictions']['clean']['ls']['rinex']['source_observation'])
    clean = observations(clean_path)
    cleanfits = {w: json.loads((OUT/f'clean_{w}.json').read_text()) for w in range(27,80)}
    output = dict(scope='post-hoc, one inspected seed; centered clean difference is a proxy, not truth',
                  zero=seal['zero'], scenes={})
    for scene in SCENES:
        obs = observations(Path(seal['predictions'][scene]['ls']['rinex']['source_observation']))
        rows = []
        for w in range(40,80):
            if w not in clean or w not in obs:
                continue
            ref, mix = clean[w], obs[w]
            dt = (mix['epoch']-ref['epoch']).total_seconds()
            shared = sorted(set(ref['sats']) & set(mix['sats']))
            # Align observed ranges approximately with measured Doppler;
            # remove common clock offset using unaffected PRNs, not G09.
            differences = {s: mix['sats'][s][0]-ref['sats'][s][0]
                           + (299792458/1575.42e6)*(mix['sats'][s][2] or 0)*dt
                           for s in shared}
            controls = [differences[s] for s in shared if s != 'G09']
            common = float(np.median(controls))
            proxy = differences['G09']-common
            fit = json.loads((OUT/f'{scene}_{w}.json').read_text())
            for method in ('ls','huber'):
                folds = fit[method]
                tau = float(np.mean([f['tau_chip'] for f in folds]))
                candidate = -(tau-seal['zero'][method])*helper.CHIP_METRES
                current_clean_tau = float(np.mean([f['tau_chip'] for f in cleanfits[w][method]]))
                row = dict(window=w, method=method, proxy_uncorrected_m=proxy,
                           desired_proxy_correction_m=-proxy, candidate_m=candidate,
                           proxy_after_candidate_m=proxy+candidate,
                           proxy_improves=abs(proxy+candidate)<abs(proxy)-1e-9,
                           fold_disagreement_m=abs(folds[0]['tau_chip']-folds[1]['tau_chip'])*helper.CHIP_METRES,
                           cv_gain_min=min(f['gain'] for f in folds),
                           boundary=any(f['boundary'] for f in folds),
                           max_rho=max(f['rho'] for f in folds),
                           timestamp_delta_s=dt, common_offset_m=common,
                           control_scatter_m=float(np.std(np.asarray(controls)-common)),
                           paired_clean_fit_correction_m=-(tau-current_clean_tau)*helper.CHIP_METRES)
                for branch in (method,method+'_cap'):
                    decisions = json.loads((OUT/scene/branch/'decisions.json').read_text())
                    decision = next(x for x in decisions if x['window']==w)
                    row[branch+'_active'] = decision['active']
                rows.append(row)
        stats = {}
        for method in ('ls','huber'):
            rr = [r for r in rows if r['method']==method]
            stats[method] = dict(proxy_uncorrected=summary([r['proxy_uncorrected_m'] for r in rr]),
                control_scatter=summary([r['control_scatter_m'] for r in rr]),
                fold_agreement_count=sum(r['fold_disagreement_m']<1e-8 for r in rr),
                candidate_improves_count=sum(r['proxy_improves'] for r in rr),
                calibration_tau=summary([np.mean([f['tau_chip'] for f in cleanfits[w][method]])*helper.CHIP_METRES for w in range(27,40)]),
                heldout_clean_tau=summary([np.mean([f['tau_chip'] for f in cleanfits[w][method]])*helper.CHIP_METRES for w in range(40,80)]),
                active_capped=[r for r in rr if r[method+'_cap_active']])
        cell = 'd0p05_a0p20_p090' if scene == 'clean' else scene
        directory = ('pvt_gated_shared_ratio_complete_hostgpu' if cell.endswith('p000')
                     else f'pvt_gated_shared_ratio_{cell}_hostgpu')
        source = BASE/directory/('clean' if scene == 'clean' else 'mixed')/'correlator_scan/epochs.csv'
        # Select one tap per epoch; no averaging or resampling of the source.
        awk = ('NR>1 && $13==0 {w=$1;if(!(w in lo))lo[w]=$3;hi[w]=$3;n[w]++;'
               'e=$4/25000-$3;if(e<0)e=-e;if(e>err)err=e}'
               'END{for(w in lo)print w,lo[w],hi[w],n[w];print "error",err+0}')
        timing = subprocess.check_output(['awk','-F,',awk,str(source)],text=True)
        windows = {}
        for line in timing.splitlines():
            f = line.split()
            if f[0] == 'error':
                time_error = float(f[1])
            else:
                w = int(f[0]); start,end,n = float(f[1])/1000,float(f[2])/1000,int(f[3])
                windows[w] = dict(first_epoch_offset_s=start-w,
                                  last_epoch_offset_s=end-w, epochs=n)
        output['scenes'][scene] = dict(stats=stats, windows=rows,
            timing=dict(sample_time_max_error_ms=time_error, windows=windows))
    (DEST/'DIAGNOSTIC.json').write_text(json.dumps(output,indent=2))
    for scene, data in output['scenes'].items():
        print(scene, json.dumps(data['stats']))


if __name__ == '__main__':
    main()
