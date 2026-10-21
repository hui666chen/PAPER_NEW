"""Only tau is refined, conditional on the saved LS separation and complex rho."""
import argparse
import copy
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar
from run import ROOT, BASE, BANK, SCENES, old, ResponseSpec, load_response_model
from fit import evaluate
from score import decisions, helper

PREVIOUS = ROOT/'results/huber_binary_gate_20261020'
OUTPUT = ROOT/'results/tau_only_local_refinement_20261021'


def refine(objective, center, halfwidth):
    lo,hi = max(-.4,center-halfwidth), min(.4,center+halfwidth)
    initial = float(objective(center))
    optimum = minimize_scalar(objective, bounds=(lo,hi), method='bounded',
                              options=dict(xatol=1e-7,maxiter=100))
    candidates = [(initial,center),(float(objective(lo)),lo),(float(objective(hi)),hi)]
    if optimum.success and np.isfinite(optimum.fun):
        candidates.append((float(optimum.fun),float(optimum.x)))
    cost,tau = min(candidates)
    return dict(tau_chip=tau,train_rss_before=initial,train_rss_after=cost,
                local_interval=[lo,hi],local_edge=bool(min(abs(tau-lo),abs(tau-hi))<1e-6),
                optimizer_success=bool(optimum.success),nfev=int(optimum.nfev))


def source(scene):
    cell = 'd0p05_a0p20_p090' if scene=='clean' else scene
    directory = 'pvt_gated_shared_ratio_complete_hostgpu' if cell.endswith('p000') else f'pvt_gated_shared_ratio_{cell}_hostgpu'
    return BASE/directory/('clean' if scene=='clean' else 'mixed')


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--scenes',nargs='+',default=SCENES)
    args=ap.parse_args();OUTPUT.mkdir(exist_ok=True)
    protocol=ROOT/'protocols/TAU_ONLY_LOCAL_REFINEMENT_20261021.md'
    manifest=dict(protocol_sha256=helper.sha256(protocol),bank_sha256=helper.sha256(BANK),
                  script_sha256=helper.sha256(Path(__file__)),truth_accessed=False)
    (OUTPUT/'MANIFEST.json').write_text(json.dumps(manifest,indent=2))
    model=load_response_model(BANK,expected_spec=ResponseSpec(prn=9))
    scale=np.linalg.norm(model.exact([0])[0])
    for scene in args.scenes:
        saved={r['window']:r['fit'] for r in json.loads((source(scene)/'gate/gate_records.json').read_text())}
        wanted=tuple(w for w in range(27,80) if not (OUTPUT/f'{scene}_{w}.json').exists())
        if not wanted: continue
        data=old.load_windows(source(scene)/'correlator_scan/epochs.csv',wanted)
        for w in wanted:
            previous=json.loads((PREVIOUS/f'{scene}_{w}.json').read_text())
            record=copy.deepcopy(previous)
            y,weights,_=data[w]
            refinements=[]
            for i,key in enumerate(('A_to_B','B_to_A')):
                train=np.arange(i,100,2); test=np.arange(1-i,100,2)
                fit=saved[w][key]
                rho=complex(fit['h1']['rho_real'],fit['h1']['rho_imag'])
                sep=fit['h1']['separation_chip']
                def template(tau,two=True):
                    response=model.exact([-tau,sep-tau] if two else [-tau])/scale
                    return response[0]+rho*response[1] if two else response[0]
                f0=refine(lambda t:evaluate(y[train],weights[train],template(t,False)),fit['h0']['tau_chip'],.005)
                f1=refine(lambda t:evaluate(y[train],weights[train],template(t)),fit['h1']['tau_chip'],.01)
                c0=evaluate(y[test],weights[test],template(f0['tau_chip'],False))
                c1=evaluate(y[test],weights[test],template(f1['tau_chip']))
                item=record['ls'][i]
                item.update(tau_chip=f1['tau_chip'],gain=(c0-c1)/max(c0,1e-30),
                            h0_cost=c0,h1_cost=c1,ordinary_rss=c1,
                            boundary=bool(np.isclose(abs(f1['tau_chip']),.4) or
                                          np.isclose(sep,.04) or np.isclose(sep,1.20)))
                refinements.append(dict(h0=f0,h1=f1,rho_real=rho.real,rho_imag=rho.imag,
                                        separation_chip=sep))
            # Huber is intentionally not part of this experiment.
            record.pop('huber',None)
            record['refinements']=refinements
            record['previous_result_sha256']=helper.sha256(PREVIOUS/f'{scene}_{w}.json')
            (OUTPUT/f'{scene}_{w}.json').write_text(json.dumps(record,indent=2))
            print(scene,w,flush=True)
    if all((OUTPUT/f'{s}_{w}.json').exists() for s in SCENES for w in range(27,80)):
        summarize()


def summarize():
    output={}
    cleanold=[json.loads((PREVIOUS/f'clean_{w}.json').read_text()) for w in range(27,80)]
    cleannew=[json.loads((OUTPUT/f'clean_{w}.json').read_text()) for w in range(27,80)]
    def tau(r):return np.mean([x['tau_chip'] for x in r['ls']])
    zero_old=float(np.median([tau(r) for r in cleanold[:13]]))
    zero_new=float(np.median([tau(r) for r in cleannew[:13]]))
    clean_stats=dict(zero_old_chip=zero_old,zero_new_chip=zero_new,
        old_std_m=float(np.std([tau(r) for r in cleanold[13:]])*helper.CHIP_METRES),
        new_std_m=float(np.std([tau(r) for r in cleannew[13:]])*helper.CHIP_METRES),
        old_centered_rms_m=float(np.sqrt(np.mean([(tau(r)-zero_old)**2 for r in cleanold[13:]]))*helper.CHIP_METRES),
        new_centered_rms_m=float(np.sqrt(np.mean([(tau(r)-zero_new)**2 for r in cleannew[13:]]))*helper.CHIP_METRES))
    for scene in SCENES:
        before=[json.loads((PREVIOUS/f'{scene}_{w}.json').read_text()) for w in range(27,80)]
        after=[json.loads((OUTPUT/f'{scene}_{w}.json').read_text()) for w in range(27,80)]
        def stats(records):
            diff=np.asarray([abs(r['ls'][0]['tau_chip']-r['ls'][1]['tau_chip'])*helper.CHIP_METRES for r in records[13:]])
            gated=decisions(records,'ls_cap',zero_old if records is before else zero_new)
            return dict(fold_disagreement_median_m=float(np.median(diff)),
                        fold_disagreement_rms_m=float(np.sqrt(np.mean(diff**2))),
                        heldout_rss_median=float(np.median([sum(x['h1_cost'] for x in r['ls']) for r in records[13:]])),
                        active=sum(r['active'] for r in gated),
                        active_evaluation=sum(r['active'] for r in gated if r['window']>=40),
                        boundary=sum(any(x['boundary'] for x in r['ls']) for r in records),
                        decisions=gated)
        output[scene]=dict(before=stats(before),after=stats(after),
            local_edge_folds=sum(f['h1']['local_edge'] for r in after for f in r['refinements']),
            optimizer_failures=sum(not f[k]['optimizer_success'] for r in after for f in r['refinements'] for k in ('h0','h1')),
            maximum_training_rss_increase=max(f[k]['train_rss_after']-f[k]['train_rss_before'] for r in after for f in r['refinements'] for k in ('h0','h1')))
    (OUTPUT/'SUMMARY.json').write_text(json.dumps(dict(clean=clean_stats,scenes=output),indent=2))
    print('SUMMARY',json.dumps(clean_stats),flush=True)


if __name__=='__main__': main()
