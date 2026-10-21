"""Post-hoc full original LS grid profile for the sole nonzero capped window."""
import json
import numpy as np
from run import old, BASE, BANK, ResponseSpec, load_response_model
from diagnose import DEST


def main():
    source = BASE/'pvt_gated_shared_ratio_d0p05_a0p50_p090_hostgpu/mixed/correlator_scan/epochs.csv'
    y, weights, _ = old.load_windows(source, (66,))[66]
    delay = np.arange(-.4,1.60001,.01)
    spec = ResponseSpec(prn=9)
    dct = load_response_model(BANK,expected_spec=spec).exact(delay).T
    rows, taus, seps = old.build_candidates(np.column_stack((dct,np.ones(401),spec.tap_grid)),delay)
    dct /= np.linalg.norm(dct[:,40])
    gram = dct.conj().T @ dct
    d,r = rows[:,0],rows[:,1]
    tau_values, inverse = np.unique(taus, return_inverse=True)
    output = {}
    for label, indices in [('even',np.arange(0,100,2)),('odd',np.arange(1,100,2)),('all',np.arange(100))]:
        corr = y[indices] @ dct.conj()
        energy = np.sum(weights[indices,None]*abs(y[indices])**2)
        profile = np.full(len(tau_values), np.inf)
        for rho in old.make_rho_grid():
            norm = gram[d,d].real + abs(rho)**2*gram[r,r].real+2*(rho*gram[d,r]).real
            rhs = corr[:,d]+rho.conjugate()*corr[:,r]
            cost = energy-np.sum(weights[indices,None]*abs(rhs)**2/norm,axis=0)
            np.minimum.at(profile,inverse,cost)
        best = int(np.argmin(profile))
        output[label] = dict(tau_chip=tau_values.tolist(),rss=profile.tolist(),
            delta_rss=(profile-profile[best]).tolist(),winner_tau=float(tau_values[best]),
            note='Original LS grid only, nuisance ratio/separation/gain profiled; not a confidence interval.')
    DEST.mkdir(exist_ok=True)
    (DEST/'WINDOW66_LS_PROFILE.json').write_text(json.dumps(output,indent=2))
    for label,item in output.items():
        print(label, 'winner',item['winner_tau'],[(round(t,2),round(j,6)) for t,j in zip(item['tau_chip'],item['delta_rss']) if .059<t<.121])


if __name__=='__main__':
    main()
