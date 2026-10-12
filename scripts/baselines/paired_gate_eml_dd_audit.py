#!/usr/bin/env python3
"""Audit paired PRN9 gated-dual-path vs EML/Double-Delta PVT results.

The script uses the already retained five-seed PRN9 batch. It checks source
RINEX/navigation hashes, verifies that the classical delays can be reproduced
from the retained 401-tap window-average scans, checks that only G09 C1C was
modified in each derived observation file, and compares seed-paired PVT RMSE.
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import sys

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[2]
WORKSPACE_ROOT = REPO_ROOT.parent
EXPERIMENT_ROOT = WORKSPACE_ROOT / "gps_l1_prn09_comparable_v3"
GATE_ROOT = EXPERIMENT_ROOT / "persistence_pvt_confirmation"
BASELINE_ROOT = EXPERIMENT_ROOT / "baseline_external_20261012/pvt_external_v2"
RESPONSE_ROOT = EXPERIMENT_ROOT / "confirmation_observation"
BASELINE_OBSERVATION_ROOT = EXPERIMENT_ROOT / "baseline_external_20261012"
OUTPUT_JSON = REPO_ROOT / "results/baselines/PRN09_GATE_VS_CLASSICAL_PAIRED_20261012.json"
OUTPUT_MD = REPO_ROOT / "results/baselines/PRN09_GATE_VS_CLASSICAL_PAIRED_20261012.md"
SEED_BOOTSTRAP = 2026101099
BOOTSTRAP_REPLICATES = 200_000
SEED_SCENES = ("clean", "amp_0p05", "amp_0p10")
CLASSICAL_METHODS = ("eml", "double_delta")
GATE_METHOD = "gated_persistence_shared_ratio"
METRICS = (
    "horizontal_rmse_m",
    "vertical_rmse_m",
    "rmse_3d_m",
    "p95_3d_m",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def mean(values: list[float]) -> float:
    return float(np.mean(np.asarray(values, dtype=float)))


def paired_interval(values: np.ndarray, resample_indices: np.ndarray) -> list[float]:
    draws = values[resample_indices]
    means = np.mean(draws, axis=1)
    return [float(x) for x in np.quantile(means, [0.025, 0.975])]


def audit_rinex_pair(reference: Path, candidate: Path, satellite: str = "G09") -> dict:
    ref_lines = reference.read_text(encoding="ascii").splitlines()
    out_lines = candidate.read_text(encoding="ascii").splitlines()
    if len(ref_lines) != len(out_lines):
        raise RuntimeError(f"RINEX line-count mismatch: {candidate}")
    changed = 0
    bad_lines: list[int] = []
    for line_number, (left, right) in enumerate(zip(ref_lines, out_lines), start=1):
        if left == right:
            continue
        changed += 1
        if (
            left[:3] != satellite
            or right[:3] != satellite
            or left[17:19] != right[17:19]
            or left[19:] != right[19:]
        ):
            bad_lines.append(line_number)
    if bad_lines:
        raise RuntimeError(
            f"{candidate}: differences outside {satellite} C1C at {bad_lines[:5]}"
        )
    return {"changed_lines": changed, "only_target_c1c_changed": True}


def main() -> None:
    # Import the exact frozen EML/DD implementations used for the original run.
    sys.path.insert(0, str(REPO_ROOT / "scripts/baselines"))
    import correlator_baselines as baselines  # noqa: PLC0415

    gate_summary = read_json(GATE_ROOT / "PVT_CONFIRMATION_SUMMARY.json")
    baseline_summary = read_json(BASELINE_ROOT / "EXTERNAL_BASELINE_PVT_SUMMARY.json")
    seeds = [int(seed) for seed in gate_summary["seeds"]]
    if seeds != [int(seed) for seed in baseline_summary["seeds"]]:
        raise RuntimeError("gate and classical baseline seed lists differ")

    gate_metrics: dict[tuple[int, str, str], dict[str, str]] = {}
    baseline_metrics: dict[tuple[int, str, str], dict[str, str]] = {}
    seed_provenance = []
    input_delay_checks = []
    rinex_checks = []
    single_solution_checks = []

    for gate_seed, baseline_seed in zip(
        gate_summary["per_seed"], baseline_summary["per_seed"]
    ):
        seed = int(gate_seed["seed"])
        if seed != int(baseline_seed["seed"]):
            raise RuntimeError("seed ordering mismatch")
        gate_rinex = gate_seed["rinex"]
        baseline_rinex = baseline_seed["rinex"]
        same_source = (
            gate_rinex["source_observation_sha256"]
            == baseline_rinex["source_observation_sha256"]
        )
        same_navigation = (
            gate_rinex["source_navigation_sha256"]
            == baseline_rinex["source_navigation_sha256"]
        )
        if not same_source or not same_navigation:
            raise RuntimeError(f"seed {seed}: source observation/navigation mismatch")

        gate_common = GATE_ROOT / f"seed_{seed}/common_observable"
        baseline_common = BASELINE_ROOT / f"seed_{seed}/common_observable"
        gate_single = gate_common / "observations/common_single.26O"
        baseline_single = baseline_common / "observations/common_single.26O"
        if sha256(gate_single) != sha256(baseline_single):
            raise RuntimeError(f"seed {seed}: common single-path RINEX differs")
        if sha256(gate_common / "navigation/common.26N") != sha256(
            baseline_common / "navigation/common.26N"
        ):
            raise RuntimeError(f"seed {seed}: copied navigation file differs")

        gate_pos = (gate_common / "solutions/single.pos").read_text().splitlines()
        baseline_pos = (baseline_common / "solutions/single.pos").read_text().splitlines()
        strip_input_path = lambda lines: [line for line in lines if not line.startswith("% inp file")]
        same_solution = strip_input_path(gate_pos) == strip_input_path(baseline_pos)
        if not same_solution:
            raise RuntimeError(f"seed {seed}: single-method PVT solution differs")
        single_solution_checks.append({"seed": seed, "same_numerical_solution": True})

        for method in (GATE_METHOD,):
            for row in read_csv(GATE_ROOT / f"seed_{seed}/metrics.csv"):
                if row["method"] == method:
                    gate_metrics[(seed, row["scene"], method)] = row
        for row in read_csv(BASELINE_ROOT / f"seed_{seed}/metrics.csv"):
            if row["method"] in CLASSICAL_METHODS:
                baseline_metrics[(seed, row["scene"], row["method"])] = row

        for method in CLASSICAL_METHODS:
            check = audit_rinex_pair(
                baseline_common / "observations/common_single.26O",
                baseline_common / f"observations/common_{method}.26O",
            )
            rinex_checks.append({"seed": seed, "method": method, **check})
        check = audit_rinex_pair(
            gate_common / "observations/common_single.26O",
            gate_common / f"observations/common_{GATE_METHOD}.26O",
        )
        rinex_checks.append({"seed": seed, "method": GATE_METHOD, **check})

        per_seed_delay_errors = []
        for scene in SEED_SCENES:
            response_file = (
                RESPONSE_ROOT / f"seed_{seed}/open_loop/{scene}/"
                "frozen_open_loop_prn9_responses.csv"
            )
            baseline_file = BASELINE_OBSERVATION_ROOT / f"seed_{seed}_{scene}.csv"
            scans = baselines.read_epoch_scans(response_file)
            estimates = {int(row["window"]): row for row in read_csv(baseline_file)}
            if set(estimates) != {scan.window for scan in scans}:
                raise RuntimeError(f"seed {seed}/{scene}: EML/DD window sets differ")
            for scan in scans:
                saved = estimates[scan.window]
                eml = float(baselines.eml_estimate(scan, baselines.BaselineConfig())["tau_chip"])
                dd = float(
                    baselines.double_delta_estimate(scan, baselines.BaselineConfig())["tau_chip"]
                )
                per_seed_delay_errors.extend([
                    abs(eml - float(saved["eml_tau_chip"])),
                    abs(dd - float(saved["dd_tau_chip"])),
                ])
        max_delay_error = float(max(per_seed_delay_errors, default=0.0))
        if max_delay_error > 1e-12:
            raise RuntimeError(f"seed {seed}: baseline delay reproduction failed")
        input_delay_checks.append({
            "seed": seed,
            "scenes": len(SEED_SCENES),
            "windows_per_scene": 16,
            "max_abs_delay_difference_chip": max_delay_error,
        })
        seed_provenance.append({
            "seed": seed,
            "same_source_observation_sha256": same_source,
            "same_source_navigation_sha256": same_navigation,
            "same_common_single_rinex": True,
            "same_common_navigation": True,
        })

    if len(gate_metrics) != len(seeds) * len(SEED_SCENES):
        raise RuntimeError("gate metrics are incomplete")
    if len(baseline_metrics) != len(seeds) * len(SEED_SCENES) * len(CLASSICAL_METHODS):
        raise RuntimeError("classical baseline metrics are incomplete")

    rng = np.random.default_rng(SEED_BOOTSTRAP)
    resample_indices = rng.integers(
        0, len(seeds), size=(BOOTSTRAP_REPLICATES, len(seeds))
    )
    results = []
    for scene in SEED_SCENES:
        method_results = {}
        for method, table in ((GATE_METHOD, gate_metrics),):
            rows = [table[(seed, scene, method)] for seed in seeds]
            method_results[method] = {
                "seed_count": len(rows),
                "mean_pvt_availability": mean([float(row["pvt_availability"]) for row in rows]),
                "mean_metrics": {
                    metric: mean([float(row[metric]) for row in rows]) for metric in METRICS
                },
                "per_seed_pvt_availability": [float(row["pvt_availability"]) for row in rows],
            }
        for method in CLASSICAL_METHODS:
            rows = [baseline_metrics[(seed, scene, method)] for seed in seeds]
            method_results[method] = {
                "seed_count": len(rows),
                "mean_pvt_availability": mean([float(row["pvt_availability"]) for row in rows]),
                "mean_metrics": {
                    metric: mean([float(row[metric]) for row in rows]) for metric in METRICS
                },
                "per_seed_pvt_availability": [float(row["pvt_availability"]) for row in rows],
            }

        paired = {}
        for method in CLASSICAL_METHODS:
            paired[method] = {}
            for metric in METRICS:
                gate_values = np.asarray([
                    float(gate_metrics[(seed, scene, GATE_METHOD)][metric]) for seed in seeds
                ])
                baseline_values = np.asarray([
                    float(baseline_metrics[(seed, scene, method)][metric]) for seed in seeds
                ])
                # Positive means the gated method has lower error than this baseline.
                improvement = baseline_values - gate_values
                paired[method][metric] = {
                    "estimand": "baseline_metric_minus_gated_metric; positive favors gated",
                    "per_seed_improvement_m": [float(value) for value in improvement],
                    "mean_improvement_m": float(np.mean(improvement)),
                    "paired_seed_bootstrap_95pct_ci_m": paired_interval(
                        improvement, resample_indices
                    ),
                    "gate_better_tie_worse": {
                        "better": int(np.sum(improvement > 0.0)),
                        "tie": int(np.sum(improvement == 0.0)),
                        "worse": int(np.sum(improvement < 0.0)),
                    },
                }
        results.append({
            "scene": scene,
            "mean_per_seed_metrics": method_results,
            "paired_gated_vs_classical": paired,
        })

    payload = {
        "status": "matched_prn09_gate_eml_double_delta_paired_audit_complete",
        "created_date": "2026-10-12",
        "target_prn": 9,
        "seeds": seeds,
        "independent_unit": "noise seed (n=5); 16 scored seconds per seed and scene",
        "scenes": list(SEED_SCENES),
        "methods": {
            GATE_METHOD: "frozen shared-ratio dual-path with 3-of-5 persistence gate",
            "eml": "frozen EML, 0.10-chip early-late spacing",
            "double_delta": "frozen normalized Double-Delta, inner 0.10 / outer 0.50 chip",
        },
        "matched_design_audit": {
            "same_seed_sets": True,
            "same_source_observation_hash_per_seed": all(
                item["same_source_observation_sha256"] for item in seed_provenance
            ),
            "same_navigation_hash_per_seed": all(
                item["same_source_navigation_sha256"] for item in seed_provenance
            ),
            "same_common_single_observation_and_navigation_files": True,
            "only_g09_c1c_changed_in_method_observation_files": True,
            "single_method_pvt_solution_numerically_identical": True,
            "eml_dd_recomputed_from_retained_401_tap_window_average": True,
            "max_eml_dd_delay_reproduction_error_chip": max(
                item["max_abs_delay_difference_chip"] for item in input_delay_checks
            ),
            "baseline_estimator_input_level": "one coherent 401-tap average per one-second window derived from the same 100 epochs",
            "same_calibration_information_path": False,
            "calibration_asymmetry": {
                "gated_dual_path": "independently calibrated response bank and frozen existence threshold; no per-seed clean-window median is subtracted from the tau correction",
                "eml_and_double_delta": "per-seed, per-method median delay from the same seed's clean windows is subtracted as a zero point",
                "implication": "the result compares the frozen operational pipelines, but does not isolate estimator choice under identical zero-point calibration information",
            },
            "seed_provenance": seed_provenance,
            "rinex_field_audit": rinex_checks,
            "single_solution_audit": single_solution_checks,
            "delay_reproduction_audit": input_delay_checks,
        },
        "uncertainty": {
            "method": "paired percentile bootstrap over the five independent seeds",
            "replicates": BOOTSTRAP_REPLICATES,
            "random_seed": SEED_BOOTSTRAP,
            "small_sample_note": "n=5; intervals are exploratory and should not be treated as broad population guarantees",
        },
        "results": results,
        "interpretation": [
            "In clean data, the gated method and both classical correlators have similar 3D RMSE; both paired 95% intervals include zero.",
            "At amplitude ratio 0.05 and 0.10, EML and Double-Delta have lower 3D PVT RMSE than the gated dual-path method in all five paired seeds.",
            "This is a matched common-observation PVT comparison of the frozen operational pipelines on one PRN9 geometry and these three simulated scenes; because zero-point calibration differs, it does not isolate estimator choice under identical calibration information and does not establish universal ranking across PRNs, geometries, or front ends.",
        ],
        "limitations": [
            "The external baselines use a coherent one-second average of 100 epoch correlator responses; the proposed gate uses its frozen cross-fit/persistence procedure on the same underlying window data. This is matched at signal, seed, geometry, and PVT level, but it is not a per-epoch discriminator implementation comparison.",
            "Only five independent noise seeds are available.",
            "The PVT score uses a fixed nine-satellite geometry with only G09 C1C corrected; it is not a field trial or multi-PRN correction experiment.",
            "Double-Delta has one unavailable PVT epoch in the amplitude-ratio 0.10 scene (mean availability 98.75%); missing solutions are not imputed.",
        ],
    }

    lines = [
        "# PRN9 paired PVT comparison: gated dual-path vs EML and Double-Delta",
        "",
        "Date: 2026-10-12  ",
        "Status: matched paired re-analysis of the retained five-seed batch; no new IQ generation or GPU run was needed.",
        "",
        "## Design audit",
        "",
        "The gate and classical branches use the same five noise seeds, the same source observation and navigation files per seed, the same nine-satellite geometry, and the same PVT solver/configuration. Each derived RINEX stream changes only G09 C1C. The underlying common single-path PVT output is numerically identical between the two runs after excluding input-path comments. Recomputing EML and Double-Delta from the retained 401-tap window averages reproduces every stored delay exactly (maximum absolute difference 0 chip).",
        "",
        "Each seed contributes 16 scored seconds per scene. The statistical unit is the independent noise seed (n=5). The paired effect is `baseline 3D RMSE − gated 3D RMSE`; positive values favor the gate. Confidence intervals are percentile bootstrap intervals over the five paired seed differences (200,000 resamples, seed 2026101099).",
        "",
        "## 3D PVT RMSE",
        "",
        "| Scene | Gated dual-path | EML | Gate improvement vs EML (EML−gate; 95% CI) | Double-Delta | Gate improvement vs DD (DD−gate; 95% CI) |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for result in results:
        scene = result["scene"]
        methods = result["mean_per_seed_metrics"]
        paired = result["paired_gated_vs_classical"]
        gate_rmse = methods[GATE_METHOD]["mean_metrics"]["rmse_3d_m"]
        eml_rmse = methods["eml"]["mean_metrics"]["rmse_3d_m"]
        dd_rmse = methods["double_delta"]["mean_metrics"]["rmse_3d_m"]
        eml_pair = paired["eml"]["rmse_3d_m"]
        dd_pair = paired["double_delta"]["rmse_3d_m"]
        lines.append(
            f"| {scene} | {gate_rmse:.3f} | {eml_rmse:.3f} | "
            f"{eml_pair['mean_improvement_m']:+.3f} "
            f"[{eml_pair['paired_seed_bootstrap_95pct_ci_m'][0]:+.3f}, "
            f"{eml_pair['paired_seed_bootstrap_95pct_ci_m'][1]:+.3f}] | "
            f"{dd_rmse:.3f} | {dd_pair['mean_improvement_m']:+.3f} "
            f"[{dd_pair['paired_seed_bootstrap_95pct_ci_m'][0]:+.3f}, "
            f"{dd_pair['paired_seed_bootstrap_95pct_ci_m'][1]:+.3f}] |"
        )
    lines += [
        "",
        "All RMSE values are metres and are means of per-seed RMSEs. For the paired-improvement columns, a negative number means the classical baseline is better.",
        "",
        "Calibration is not fully symmetric: the classical EML/Double-Delta branches subtract a per-seed, per-method zero point estimated as the median delay over that same seed's clean windows. The gated dual-path branch uses its independently calibrated response bank and frozen existence threshold, without subtracting a per-seed clean-window median. Thus this table compares the frozen end-to-end pipelines, but it is not a pure estimator-only comparison with identical zero-point calibration information.",
        "",
        "## Availability and interpretation",
        "",
        "| Scene | Gate availability | EML availability | Double-Delta availability | Gate paired record vs EML / DD (better–tie–worse) |",
        "|---|---:|---:|---:|---:|",
    ]
    for result in results:
        scene = result["scene"]
        methods = result["mean_per_seed_metrics"]
        paired = result["paired_gated_vs_classical"]
        counts = []
        for method in CLASSICAL_METHODS:
            stat = paired[method]["rmse_3d_m"]["gate_better_tie_worse"]
            counts.append(f"{stat['better']}–{stat['tie']}–{stat['worse']}")
        lines.append(
            f"| {scene} | {methods[GATE_METHOD]['mean_pvt_availability']:.3f} | "
            f"{methods['eml']['mean_pvt_availability']:.3f} | "
            f"{methods['double_delta']['mean_pvt_availability']:.3f} | "
            f"{counts[0]} / {counts[1]} |"
        )
    lines += [
        "",
        "The result is clear for this test condition: the gate is comparable to the classical correlators in clean data, but EML and Double-Delta substantially reduce 3D PVT RMSE more in both multipath scenes. At amplitude ratio 0.05, the gate is worse by about 1.86 m; at 0.10, it is worse by about 2.79–2.89 m. All five paired seeds favor each classical baseline in both multipath scenes. This does not show that the GPU software receiver itself is worse; it shows that the current PRN9 gated correction has not beaten these classical multipath controls.",
        "",
        "The comparison is matched at the IQ-derived observation, seed, geometry, and PVT levels. It remains a window-level estimator comparison: EML/Double-Delta use a coherent one-second 401-tap average, whereas the gate uses its frozen cross-fit/persistence decision on the same underlying window data. Only five seeds and one PRN9 geometry are covered, so do not generalize beyond this setup. Because calibration information is asymmetric, treat the result as a serious operational warning and not the final estimator-only ranking; a strict method-only claim needs a follow-up with one independent clean calibration set and the same zero-point convention for every method, without retuning the frozen gate.",
        "",
        "Machine-readable values and audit details: `PRN09_GATE_VS_CLASSICAL_PAIRED_20261012.json`.",
        "",
    ]

    OUTPUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_JSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    OUTPUT_MD.write_text("\n".join(lines), encoding="utf-8")
    print(f"JSON={OUTPUT_JSON}")
    print(f"MARKDOWN={OUTPUT_MD}")
    for result in results:
        print(json.dumps({
            "scene": result["scene"],
            "gate_3d_rmse": result["mean_per_seed_metrics"][GATE_METHOD]["mean_metrics"]["rmse_3d_m"],
            "eml_3d_rmse": result["mean_per_seed_metrics"]["eml"]["mean_metrics"]["rmse_3d_m"],
            "double_delta_3d_rmse": result["mean_per_seed_metrics"]["double_delta"]["mean_metrics"]["rmse_3d_m"],
            "gate_minus_eml_paired_improvement": result["paired_gated_vs_classical"]["eml"]["rmse_3d_m"],
            "gate_minus_dd_paired_improvement": result["paired_gated_vs_classical"]["double_delta"]["rmse_3d_m"],
        }, sort_keys=True))


if __name__ == "__main__":
    main()
