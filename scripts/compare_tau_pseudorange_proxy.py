"""Compare direct-path tau error against single-component baselines.

This replays the exact old and confirmatory Monte Carlo streams.  It reports
tau error in chips and metres; no absolute pseudorange/geometric truth exists
in these synthetic MC files, so metre-valued tau error is a proxy rather than
a complete receiver pseudorange error.
"""

import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import monte_carlo_resolution_baseline as mc
import run_confirmatory_outlier_mc as confirm_mc
import run_preregistered_outlier_mc as old_mc
from blind_response_experiment import write_new
from run_blind_transfer import ROOT


HERE = Path(__file__).resolve().parent
OUT = ROOT / "tau_pseudorange_proxy_comparison"
OLD_PREDICTION = ROOT / "outlier_theory_prediction_27_v2/predictions.json"
OLD_BASE = ROOT / "outlier_theory_mc_base200"
OLD_SUPPLEMENT = ROOT / "outlier_theory_mc_supplement800"
NEW_PREDICTION = (
    ROOT / "outlier_theory_confirmatory_prediction_27_gaussian_frozen/predictions.json"
)
NEW_MC = ROOT / "outlier_theory_confirmatory_mc"
SPEED_OF_LIGHT = 299_792_458.0
GPS_CA_CODE_RATE = 1_023_000.0
METRES_PER_CHIP = SPEED_OF_LIGHT / GPS_CA_CODE_RATE
OUTLIER_CHIP = .05
GATE_THRESHOLDS = (0.0, .01, .02, .05, .10, .20, .50, 1.0)


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def cell_key(delay, amp, epochs):
    return f"d{delay:.2f}_a{amp:.2f}_n{epochs:d}"


def build_single_path_models(delay_grid, dictionary):
    q = dictionary.shape[1]
    direct = np.rint((-mc.TAU - delay_grid[0]) / .01).astype(int)
    nuisance = np.column_stack((
        np.full(len(mc.TAU), q - 2), np.full(len(mc.TAU), q - 1)
    ))
    indices = np.column_stack((direct, nuisance))
    design = np.stack([dictionary[:, row] for row in indices])
    gram = np.einsum("nmc,nmd->ncd", design.conj(), design)
    inverse = np.linalg.inv(gram)
    direct_norm2 = np.sum(abs(dictionary[:, direct]) ** 2, axis=0)
    return direct, indices, inverse, direct_norm2


def cell_moment(cell, base, delay_grid, dictionary):
    delay, ratio, _ = cell
    direct_index = int(round((0 - delay_grid[0]) / .01))
    reflected_index = int(round((delay - delay_grid[0]) / .01))
    rho = ratio * np.exp(1j * mc.PHASE_RAD)
    h = dictionary[:, direct_index] + rho * dictionary[:, reflected_index]
    hcorr = dictionary.conj().T @ h
    henergy = float(np.vdot(h, h).real)
    return mc.moments(
        base, hcorr, henergy, [(1, direct_index), (rho, reflected_index)]
    )[0]


def baseline_tau(moment, direct, indices, inverse, direct_norm2):
    sub = moment[:, indices[:, :, None], indices[:, None, :]]
    explained = np.einsum("cij,rcij->rc", inverse, sub).real
    tau_single_nls = mc.TAU[np.argmax(explained, axis=1)]
    peak_energy = moment[:, direct, direct].real / direct_norm2[None, :]
    tau_energy_peak = mc.TAU[np.argmax(peak_energy, axis=1)]
    return tau_single_nls, tau_energy_peak


def append_replay(storage, cell, base, delay_grid, dictionary,
                  direct, single_indices, single_inverse, direct_norm2):
    moment = cell_moment(cell, base, delay_grid, dictionary)
    single, peak = baseline_tau(
        moment, direct, single_indices, single_inverse, direct_norm2
    )
    storage[cell]["single_nls"].append(single)
    storage[cell]["energy_peak"].append(peak)


def replay_new(delay_grid, dictionary, eigenvalues, eigenvectors, gain_pool,
               direct, single_indices, single_inverse, direct_norm2):
    prediction = json.loads(NEW_PREDICTION.read_text())["results"]
    storage = {}
    for number, row in enumerate(prediction, start=1):
        cell = (float(row["delay"]), float(row["amp"]), int(row["N_epoch"]))
        target = 1000 if float(row["P_global_pred"]) < .02 else 200
        storage[cell] = {"single_nls": [], "energy_peak": []}
        for first in range(0, target, confirm_mc.BATCH_REPEATS):
            trials = np.arange(first, min(first + confirm_mc.BATCH_REPEATS, target))
            noise_corr, noise_energy, weights = confirm_mc.synthesize_trials(
                trials, *cell, eigenvalues, eigenvectors, dictionary, gain_pool
            )
            base = mc.noise_moments(noise_corr, noise_energy, weights, cell[2])
            append_replay(
                storage, cell, base, delay_grid, dictionary, direct,
                single_indices, single_inverse, direct_norm2,
            )
        name = cell_key(*cell) + ".npz"
        with np.load(NEW_MC / name) as saved:
            two = np.asarray(saved["winner_tau"], float)
            trial_index = np.asarray(saved["trial_index"], int)
        yield dict(
            batch="new_confirmatory_27", cell=cell,
            predicted_global=float(row["P_global_pred"]),
            cost_gap_tier=row["cost_gap_tier"],
            trial_index=trial_index, two_path=two,
            single_nls=np.concatenate(storage[cell]["single_nls"]),
            energy_peak=np.concatenate(storage[cell]["energy_peak"]),
        )
        print(f"new {number}/27 {cell_key(*cell)}", flush=True)


def replay_old(delay_grid, dictionary, eigenvalues, eigenvectors, gain_pool,
               direct, single_indices, single_inverse, direct_norm2):
    prediction = json.loads(OLD_PREDICTION.read_text())["results"]
    cells = [(float(r["delay"]), float(r["amp"]), int(r["N_epoch"])) for r in prediction]
    storage = {cell: {"single_nls": [], "energy_peak": []} for cell in cells}

    rng = np.random.default_rng(old_mc.BASE_SEED)
    noise_corr, noise_energy, weights = old_mc.synthesize(
        rng, old_mc.BASE_REPEATS, eigenvalues, eigenvectors, dictionary, gain_pool
    )
    bases = {
        count: mc.noise_moments(noise_corr, noise_energy, weights, count)
        for count in sorted({cell[2] for cell in cells})
    }
    for cell in cells:
        append_replay(
            storage, cell, bases[cell[2]], delay_grid, dictionary, direct,
            single_indices, single_inverse, direct_norm2,
        )

    supplemented = {
        (float(r["delay"]), float(r["amp"]), int(r["N_epoch"]))
        for r in prediction if float(r["P_global_pred"]) < .02
    }
    rng = np.random.default_rng(old_mc.SUPPLEMENT_SEED)
    for first in range(0, old_mc.SUPPLEMENT_REPEATS, old_mc.BATCH_REPEATS):
        repeats = min(old_mc.BATCH_REPEATS, old_mc.SUPPLEMENT_REPEATS - first)
        noise_corr, noise_energy, weights = old_mc.synthesize(
            rng, repeats, eigenvalues, eigenvectors, dictionary, gain_pool
        )
        bases = {
            count: mc.noise_moments(noise_corr, noise_energy, weights, count)
            for count in sorted({cell[2] for cell in supplemented})
        }
        for cell in supplemented:
            append_replay(
                storage, cell, bases[cell[2]], delay_grid, dictionary, direct,
                single_indices, single_inverse, direct_norm2,
            )

    prediction_by_cell = {
        (float(r["delay"]), float(r["amp"]), int(r["N_epoch"])): r
        for r in prediction
    }
    for number, cell in enumerate(cells, start=1):
        name = cell_key(*cell) + ".npz"
        with np.load(OLD_BASE / name) as base_saved:
            two_parts = [np.asarray(base_saved["winner_tau"], float)]
        if cell in supplemented:
            with np.load(OLD_SUPPLEMENT / name) as supplement_saved:
                two_parts.append(np.asarray(supplement_saved["winner_tau"], float))
        two = np.concatenate(two_parts)
        row = prediction_by_cell[cell]
        yield dict(
            batch="old_27", cell=cell,
            predicted_global=float(row["P_global_pred"]),
            cost_gap_tier=row["cost_gap_tier"],
            trial_index=np.arange(len(two)), two_path=two,
            single_nls=np.concatenate(storage[cell]["single_nls"]),
            energy_peak=np.concatenate(storage[cell]["energy_peak"]),
        )
        print(f"old {number}/27 {cell_key(*cell)}", flush=True)


def metrics(values):
    values = np.asarray(values, float)
    absolute = abs(values)
    return dict(
        bias_chip=float(values.mean()),
        rmse_chip=float(np.sqrt(np.mean(values ** 2))),
        mae_chip=float(absolute.mean()),
        p95_abs_chip=float(np.quantile(absolute, .95)),
        outlier_gt_0p05=float(np.mean(absolute > OUTLIER_CHIP + 1e-12)),
        bias_m=float(values.mean() * METRES_PER_CHIP),
        rmse_m=float(np.sqrt(np.mean(values ** 2)) * METRES_PER_CHIP),
        mae_m=float(absolute.mean() * METRES_PER_CHIP),
        p95_abs_m=float(np.quantile(absolute, .95) * METRES_PER_CHIP),
    )


def cell_summary(item):
    output = dict(
        batch=item["batch"], delay=item["cell"][0], amp=item["cell"][1],
        N_epoch=item["cell"][2], repeats=len(item["two_path"]),
        predicted_global=item["predicted_global"],
        cost_gap_tier=item["cost_gap_tier"], truth_tau_chip=0.0,
    )
    for method in ("two_path", "single_nls", "energy_peak"):
        result = metrics(item[method])
        for name, value in result.items():
            output[f"{method}_{name}"] = value
    two_abs = abs(item["two_path"])
    single_abs = abs(item["single_nls"])
    peak_abs = abs(item["energy_peak"])
    output.update(
        two_vs_single_strict_win=float(np.mean(two_abs < single_abs)),
        two_vs_single_tie=float(np.mean(two_abs == single_abs)),
        two_vs_peak_strict_win=float(np.mean(two_abs < peak_abs)),
        two_vs_peak_tie=float(np.mean(two_abs == peak_abs)),
        two_vs_single_rmse_ratio=(
            output["two_path_rmse_chip"] / output["single_nls_rmse_chip"]
            if output["single_nls_rmse_chip"] > 0 else None
        ),
        two_vs_peak_rmse_ratio=(
            output["two_path_rmse_chip"] / output["energy_peak_rmse_chip"]
            if output["energy_peak_rmse_chip"] > 0 else None
        ),
    )
    return output


def aggregate(cell_rows, label):
    output = dict(label=label, cells=len(cell_rows))
    for method in ("two_path", "single_nls", "energy_peak"):
        rmses = np.asarray([row[f"{method}_rmse_m"] for row in cell_rows])
        maes = np.asarray([row[f"{method}_mae_m"] for row in cell_rows])
        output[f"{method}_mean_cell_rmse_m"] = float(rmses.mean())
        output[f"{method}_median_cell_rmse_m"] = float(np.median(rmses))
        output[f"{method}_mean_cell_mae_m"] = float(maes.mean())
    output["two_path_cells_better_than_single"] = int(sum(
        row["two_path_rmse_m"] < row["single_nls_rmse_m"] for row in cell_rows
    ))
    output["two_path_cells_better_than_peak"] = int(sum(
        row["two_path_rmse_m"] < row["energy_peak_rmse_m"] for row in cell_rows
    ))
    return output


def gate_curve(items):
    rows = []
    for threshold in GATE_THRESHOLDS:
        cell_rmses = []
        selected_two = 0
        for item in items:
            use_two = item["predicted_global"] <= threshold
            values = item["two_path"] if use_two else item["single_nls"]
            cell_rmses.append(metrics(values)["rmse_m"])
            selected_two += int(use_two)
        rows.append(dict(
            threshold=threshold, two_path_cells=selected_two,
            single_path_cells=len(items) - selected_two,
            mean_cell_rmse_m=float(np.mean(cell_rmses)),
            median_cell_rmse_m=float(np.median(cell_rmses)),
        ))
    return rows


def main():
    OUT.mkdir(exist_ok=False)
    (_, eigenvalues, eigenvectors, delay_grid, dictionary,
     _, _, gain_pool) = old_mc.setup()
    direct, single_indices, single_inverse, direct_norm2 = build_single_path_models(
        delay_grid, dictionary
    )
    items = list(replay_old(
        delay_grid, dictionary, eigenvalues, eigenvectors, gain_pool,
        direct, single_indices, single_inverse, direct_norm2,
    ))
    items.extend(replay_new(
        delay_grid, dictionary, eigenvalues, eigenvectors, gain_pool,
        direct, single_indices, single_inverse, direct_norm2,
    ))
    for item in items:
        n = len(item["two_path"])
        assert len(item["single_nls"]) == len(item["energy_peak"]) == n

    cell_rows = [cell_summary(item) for item in items]
    fields = list(cell_rows[0])
    with (OUT / "tau_cell_comparison.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(cell_rows)

    aggregates = []
    for batch in ("old_27", "new_confirmatory_27"):
        aggregates.append(aggregate(
            [row for row in cell_rows if row["batch"] == batch], batch
        ))
    aggregates.append(aggregate(cell_rows, "all_54_equal_cell_weight"))
    gate = gate_curve(items)
    summary = dict(
        schema=1, truth_accessed=True,
        truth_scope=("synthetic MC direct-path tau truth is zero; sealed B truth was inspected "
                     "only to determine that it contains relative path truth but no absolute "
                     "geometric pseudorange"),
        absolute_pseudorange_truth_available=False,
        reported_quantity=("direct-path tau error converted to metres; proxy for pseudorange "
                           "correction error, not full receiver pseudorange/PVT error"),
        metres_per_chip=METRES_PER_CHIP,
        estimator_grid_chip=.01,
        estimator_grid_m=.01 * METRES_PER_CHIP,
        truth_tau_chip=0.0,
        baselines=dict(
            single_nls=("one response column plus the same complex constant and linear nuisance "
                        "columns; global tau-grid least squares"),
            energy_peak=("noncoherent multi-epoch single-component normalized correlation-energy "
                         "peak; proxy only, not a recorded native receiver DLL output"),
        ),
        source_hashes=dict(
            old_prediction=sha256(OLD_PREDICTION),
            new_prediction=sha256(NEW_PREDICTION),
            old_base_manifest=sha256(OLD_BASE / "manifest.json"),
            old_supplement_manifest=sha256(OLD_SUPPLEMENT / "manifest.json"),
            new_mc_manifest=sha256(NEW_MC / "manifest.json"),
        ),
        aggregates=aggregates,
        risk_gate_curve=gate,
        cells=cell_rows,
    )
    write_new(OUT / "summary.json", summary)

    with (OUT / "aggregate_summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(aggregates[0]))
        writer.writeheader()
        writer.writerows(aggregates)
    with (OUT / "risk_gate_curve.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(gate[0]))
        writer.writeheader()
        writer.writerows(gate)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
    x = np.arange(len(cell_rows))
    axes[0].plot(x, [r["single_nls_rmse_m"] for r in cell_rows], "o", ms=4,
                 label="Single-path NLS")
    axes[0].plot(x, [r["two_path_rmse_m"] for r in cell_rows], "o", ms=4,
                 label="Two-path M0")
    axes[0].axvline(26.5, color="0.5", ls="--", lw=1)
    axes[0].set_xlabel("Cell (old 27 | new 27)")
    axes[0].set_ylabel("Direct-path tau RMSE (m)")
    axes[0].grid(alpha=.2)
    axes[0].legend()
    axes[1].plot(
        [r["threshold"] for r in gate], [r["mean_cell_rmse_m"] for r in gate], "o-"
    )
    axes[1].set_xscale("symlog", linthresh=.01)
    axes[1].set_xlabel("Use two-path if predicted global risk <= threshold")
    axes[1].set_ylabel("Equal-cell mean tau RMSE (m)")
    axes[1].grid(alpha=.2)
    fig.tight_layout()
    fig.savefig(OUT / "tau_rmse_and_risk_gate.png", dpi=180)
    fig.savefig(OUT / "tau_rmse_and_risk_gate.pdf")
    plt.close(fig)

    print(json.dumps({"aggregates": aggregates, "risk_gate_curve": gate}, indent=2))
    print(f"COMPLETE={OUT}")


if __name__ == "__main__":
    main()
