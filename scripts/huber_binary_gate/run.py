#!/usr/bin/env python3
import argparse
import hashlib
import json
import sys
import time
from pathlib import Path
import numpy as np
from fit import Grid, gain_residual, evaluate

ROOT = Path(__file__).resolve().parents[2]
PAPER = ROOT.parent
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER/"gps_l1_prn09_shared_ratio_v2"))
import shared_ratio_existence as old
from receiver_chain_response_control import ResponseSpec, load_response_model

BASE = ROOT/"results/prn09_phase_amp_delay_grid_20261019/seed_2026101922"
BANK = ROOT/"results/three_ray_25msps_pilot_v2_20261018/calibration/response_banks/bank_12.npz"
SCENES = ["clean", "d0p05_a0p20_p000", "d0p05_a0p20_p090",
          "d0p05_a0p20_p180", "d0p05_a0p50_p090"]


def describe(index, rows, taus, seps, ratios):
    k, row = divmod(index, len(rows))
    tau, sep, rho = float(taus[row]), float(seps[row]), float(abs(ratios[k]))
    return dict(tau_chip=tau, separation_chip=sep, rho=rho,
                boundary=bool(np.isclose(abs(tau), .4) or np.isclose(sep, .04)
                              or np.isclose(sep, 1.20)))


def analyze(y, weights, h0, h1, rows, taus, seps, ratios):
    result = {"ls": [], "huber": []}
    for train, test in ((np.arange(0, 100, 2), np.arange(1, 100, 2)),
                        (np.arange(1, 100, 2), np.arange(0, 100, 2))):
        initial = h0.fit(y[train], weights[train])
        r = initial["residual"]
        components = np.r_[r.real.ravel(), r.imag.ravel()]
        delta = max(1e-12, 1.345*np.median(abs(components-np.median(components)))/.67448975)
        for method in result:
            cutoff = delta if method == "huber" else None
            f0 = initial if cutoff is None else h0.fit(y[train], weights[train], cutoff)
            f1 = h1.fit(y[train], weights[train], cutoff)
            cost0 = evaluate(y[test], weights[test], h0.template(f0["index"]), cutoff)
            cost1 = evaluate(y[test], weights[test], h1.template(f1["index"]), cutoff)
            item = describe(f1["index"], rows, taus, seps, ratios)
            item.update(delta=delta, gain=(cost0-cost1)/max(cost0, 1e-30),
                        h0_cost=cost0, h1_cost=cost1,
                        ordinary_rss=evaluate(y[test], weights[test], h1.template(f1["index"])),
                        converged=f0["converged"] and f1["converged"],
                        objective=f1["objective"])
            result[method].append(item)
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scenes", nargs="+", default=SCENES)
    ap.add_argument("--windows", nargs="+", type=int, default=list(range(27, 80)))
    ap.add_argument("--out", type=Path, default=ROOT/"results/huber_binary_gate_20261020")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    spec = ResponseSpec(prn=9)
    delay = np.arange(-.4, 1.60001, .01)
    dictionary = load_response_model(BANK, expected_spec=spec).exact(delay).T
    augmented = np.column_stack((dictionary, np.ones(401), spec.tap_grid))
    rows, taus, seps = old.build_candidates(augmented, delay)
    ratios = old.make_rho_grid()
    dictionary = dictionary/np.linalg.norm(dictionary[:, 40])
    h1 = Grid(dictionary, rows[:, 0], rows[:, 1], ratios)
    indices = old.direct_indices(delay)
    h0 = Grid(dictionary, indices, indices, [0])
    manifest = dict(bank_sha256=hashlib.sha256(BANK.read_bytes()).hexdigest(),
                    windows=args.windows, scenes=args.scenes, truth_accessed=False,
                    protocol=str(ROOT/"protocols/HUBER_BINARY_GATE_EXPLORATORY_20261020.md"))
    (args.out/("manifest_"+"_".join(args.scenes)+".json")).write_text(json.dumps(manifest, indent=2))
    for scene in args.scenes:
        cell = "d0p05_a0p20_p090" if scene == "clean" else scene
        directory = ("pvt_gated_shared_ratio_complete_hostgpu" if cell == "d0p05_a0p20_p000"
                     else f"pvt_gated_shared_ratio_{cell}_hostgpu")
        source = BASE/directory/(
            "clean" if scene == "clean" else "mixed")/"correlator_scan/epochs.csv"
        pending = [w for w in args.windows if not (args.out/f"{scene}_{w}.json").exists()]
        if not pending:
            continue
        data = old.load_windows(source, tuple(pending))
        for w in pending:
            started = time.time()
            y, weights, _ = data[w]
            record = analyze(y, weights, h0, h1, rows, taus, seps, ratios)
            record.update(scene=scene, window=w, elapsed_seconds=time.time()-started)
            (args.out/f"{scene}_{w}.json").write_text(json.dumps(record, indent=2))
            print(scene, w, round(record["elapsed_seconds"], 2),
                  {m: [(round(x["rho"], 2), x["boundary"], x["converged"])
                       for x in record[m]] for m in ("ls", "huber")}, flush=True)


if __name__ == "__main__":
    main()
