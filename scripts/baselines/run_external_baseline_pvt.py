#!/usr/bin/env python3
"""Run external-baseline PVT from the retained common PRN9 observations.

The source RINEX for each seed is the already-generated zero-correction
receiver stream used by the formal PRN9 confirmation.  Only the target G09
C1C field is changed in each derived RINEX; all other satellites, fields,
epochs, and navigation data are kept byte-identical.  Baseline delays are
converted to schedules by the separate observation-only summary script.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil
import sys

import numpy as np


PAPER_ROOT = Path("/home/zoxoai/Documents/ChatGPT/paper")
PRN9_ROOT = PAPER_ROOT / "gps_l1_prn09_comparable_v3"
PVT_ROOT = PRN9_ROOT / "persistence_pvt_confirmation"
BASELINE_ROOT = PRN9_ROOT / "baseline_external_20261012"
OUTPUT_ROOT = BASELINE_ROOT / "pvt_external"

sys.path.insert(0, str(PAPER_ROOT / "gps_l1_prn30_shared_ratio_v2"))
import run_prn30_v2_pvt as base  # noqa: E402


TARGET_PRN = 9
TARGET_SATELLITE = "G09"
METHODS = ("single", "eml", "double_delta", "medll_ideal", "medll_chain")
SCENES = {
    "clean": (tuple(range(61, 77)), 68),
    "amp_0p05": (tuple(range(100, 116)), 108),
    "amp_0p10": (tuple(range(121, 137)), 129),
}
SEEDS = (2026101110, 2026101111, 2026101112, 2026101113, 2026101114)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_schedule(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def combine_schedule(
    paths: list[Path], output: Path, *, zero: bool = False
) -> dict[str, object]:
    rows_by_second: dict[int, dict[str, str]] = {}
    for path in paths:
        for row in read_schedule(path):
            second = int(row["start_sample"]) // base.SAMPLE_RATE_HZ
            if second in rows_by_second:
                raise ValueError(f"duplicate schedule second {second}: {path}")
            row = dict(row)
            row["prn"] = str(TARGET_PRN)
            if zero:
                row["correction_m"] = "0.000000000000"
            rows_by_second[second] = row
    rows = [rows_by_second[key] for key in sorted(rows_by_second)]
    if not rows:
        raise ValueError(f"empty schedule set for {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=("start_sample", "end_sample", "prn", "correction_m"),
        )
        writer.writeheader()
        writer.writerows({key: row[key] for key in writer.fieldnames} for row in rows)
    return {"path": str(output), "sha256": sha256(output), "rows": len(rows)}


def configure_base() -> None:
    base.TARGET_PRN = TARGET_PRN
    base.TARGET_SATELLITE = TARGET_SATELLITE
    base.METHODS = METHODS
    base.SCENES = SCENES
    base.EXPECTED_SATELLITES = (
        "G03", "G04", "G07", "G08", "G09", "G14", "G16", "G27", "G30"
    )
    base.EXPECTED_SATELLITE_COUNT = len(base.EXPECTED_SATELLITES)
    base.DYNAMIC_EXPECTED_SATELLITES = False


def find_single(directory: Path, suffix: str) -> Path:
    paths = sorted(directory.glob(f"*.??{suffix}"))
    if len(paths) != 1:
        raise RuntimeError(f"expected one *.{suffix} in {directory}, found {paths}")
    return paths[0]


def make_method_schedules(seed: int, run_dir: Path) -> tuple[dict[str, Path], dict]:
    source = BASELINE_ROOT / "schedules"
    schedule_paths: dict[str, Path] = {}
    audits: dict[str, object] = {}
    scene_paths = {
        scene: [
            source / f"seed_{seed}_{scene}_{method}.csv"
            for method in ("eml", "double_delta", "medll_ideal", "medll_chain")
        ]
        for scene in SCENES
    }
    for scene, paths in scene_paths.items():
        for path in paths:
            if not path.is_file():
                raise FileNotFoundError(path)
    first_scene = list(SCENES)[0]
    for method in METHODS:
        destination = run_dir / f"{method}_corrections.csv"
        if method == "single":
            all_paths = [scene_paths[first_scene][0]]
            for scene in list(SCENES)[1:]:
                all_paths.append(scene_paths[scene][0])
            audits[method] = combine_schedule(all_paths, destination, zero=True)
        else:
            all_paths = [
                source / f"seed_{seed}_{scene}_{method}.csv" for scene in SCENES
            ]
            audits[method] = combine_schedule(all_paths, destination)
        schedule_paths[method] = destination
    return schedule_paths, audits


def process_seed(seed: int, output_root: Path) -> dict:
    source_receiver = PVT_ROOT / f"seed_{seed}" / "source_receiver"
    source_obs = find_single(source_receiver, "O")
    source_nav = find_single(source_receiver, "N")
    run_dir = output_root / f"seed_{seed}"
    if run_dir.exists():
        raise RuntimeError(f"refusing to overwrite {run_dir}")
    run_dir.mkdir(parents=True)
    schedules, schedule_audit = make_method_schedules(seed, run_dir)
    rinex_audit = base.build_common_rinex(run_dir, source_obs, source_nav, schedules)
    base.solve_pvt(run_dir)
    result = score_external_seed(
        run_dir, seed,
        {"paths": {key: str(value) for key, value in schedules.items()},
         "sha256": {key: sha256(value) for key, value in schedules.items()},
         "audits": schedule_audit},
        {"source_observation_sha256": sha256(source_obs),
         "source_navigation_sha256": sha256(source_nav),
         "truth_accessed": False},
        rinex_audit,
    )
    result["truth_used_for_estimation"] = False
    result["truth_used_for_scoring_only"] = True
    return result


def score_external_seed(
    run_dir: Path, seed: int, schedules: dict, source: dict, rinex: dict
) -> dict:
    """Score every method while retaining partial PVT availability."""
    common = run_dir / "common_observable"
    solutions = {
        method: base.parse_solution(common / "solutions" / f"{method}.pos")
        for method in METHODS
    }
    metrics: list[dict[str, object]] = []
    epoch_rows: list[dict[str, object]] = []
    missing_solution_seconds: dict[str, list[int]] = {}
    for scene, (seconds, _) in SCENES.items():
        for method in METHODS:
            selected = []
            missing = []
            for second in seconds:
                row = solutions[method].get(second)
                if row is None or row["satellites"] != base.EXPECTED_SATELLITE_COUNT:
                    missing.append(int(second))
                    continue
                selected.append(row)
                epoch_rows.append({"seed": seed, "scene": scene, "method": method, **row})
            missing_solution_seconds[f"{scene}/{method}"] = missing
            if selected:
                row = {"seed": seed, **base.metric_row(scene, method, selected)}
            else:
                row = {
                    "seed": seed, "scene": scene, "method": method,
                    "epochs": 0, "satellites_min": None, "satellites_max": None,
                    "horizontal_rmse_m": None, "vertical_rmse_m": None,
                    "rmse_3d_m": None, "p95_3d_m": None,
                    "gross_error_rate": None, "pvt_availability": 0.0,
                    "mean_east_error_m": None, "mean_north_error_m": None,
                    "mean_up_error_m": None,
                }
            metrics.append(row)
    with (run_dir / "metrics.csv").open("w", newline="", encoding="utf-8") as stream:
        fields = list(metrics[0])
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(metrics)
    if epoch_rows:
        with (run_dir / "epoch_errors.csv").open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(epoch_rows[0]))
            writer.writeheader()
            writer.writerows(epoch_rows)
    return {
        "seed": seed,
        "status": "valid_with_partial_pvt_availability",
        "target_prn": TARGET_PRN,
        "invariants": {
            "one_common_observation_stream": True,
            "only_target_c1c_differs": True,
            "pvt_availability_scored_per_method": True,
            "expected_satellites": base.EXPECTED_SATELLITE_COUNT,
        },
        "missing_solution_seconds": missing_solution_seconds,
        "source_receiver": source,
        "rinex": rinex,
        "schedules": schedules,
        "metrics": metrics,
    }


def aggregate_external(results: list[dict]) -> dict:
    rows = [row for result in results for row in result["metrics"]]
    by_key = {(int(row["seed"]), row["scene"], row["method"]): row for row in rows}
    output = {"per_method_scene": {}, "paired_vs_single": {}}
    for scene in SCENES:
        for method in METHODS:
            selected = [by_key[(seed, scene, method)] for seed in SEEDS if (seed, scene, method) in by_key]
            available = [row for row in selected if row["rmse_3d_m"] is not None]
            def mean(key: str):
                values = [row[key] for row in available if row[key] is not None]
                return float(np.mean(values)) if values else None
            output["per_method_scene"][f"{scene}/{method}"] = {
                "seeds": len(selected),
                "available_seed_count": len(available),
                "horizontal_rmse_mean_m": mean("horizontal_rmse_m"),
                "vertical_rmse_mean_m": mean("vertical_rmse_m"),
                "rmse_3d_mean_m": mean("rmse_3d_m"),
                "p95_3d_mean_m": mean("p95_3d_m"),
                "gross_error_rate_mean": mean("gross_error_rate"),
                "pvt_availability_mean": float(np.mean([
                    row["pvt_availability"] for row in selected
                ])) if selected else 0.0,
            }
        for method in METHODS[1:]:
            differences = []
            for seed in SEEDS:
                candidate = by_key[(seed, scene, method)]
                single = by_key[(seed, scene, "single")]
                if candidate["rmse_3d_m"] is not None and single["rmse_3d_m"] is not None:
                    differences.append(single["rmse_3d_m"] - candidate["rmse_3d_m"])
            values = np.asarray(differences, dtype=float)
            output["paired_vs_single"][f"{scene}/{method}"] = {
                "paired_seed_count": int(values.size),
                "improvement_3d_rmse_m_mean": float(np.mean(values)) if values.size else None,
                "wins": int(np.sum(values > 0.0)),
                "ties": int(np.sum(values == 0.0)),
                "losses": int(np.sum(values < 0.0)),
                "improvement_values_m": values.tolist(),
            }
    return output


def main() -> None:
    global BASELINE_ROOT, PVT_ROOT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-root", type=Path, default=BASELINE_ROOT)
    parser.add_argument("--pvt-root", type=Path, default=PVT_ROOT)
    parser.add_argument("--output-root", type=Path, default=OUTPUT_ROOT)
    parser.add_argument("--seeds", nargs="+", type=int, default=list(SEEDS))
    args = parser.parse_args()
    BASELINE_ROOT = args.baseline_root
    PVT_ROOT = args.pvt_root
    configure_base()
    if args.output_root.exists() and any(args.output_root.iterdir()):
        raise RuntimeError(f"refusing non-empty output root: {args.output_root}")
    args.output_root.mkdir(parents=True, exist_ok=True)
    results = []
    for seed in args.seeds:
        print(f"seed={seed} stage=common_rinex_and_pvt", flush=True)
        results.append(process_seed(seed, args.output_root))
    base.CONFIRMATION_SEEDS = tuple(args.seeds)
    aggregate = aggregate_external(results)
    payload = {
        "status": "external_baseline_common_observable_pvt_complete",
        "target_prn": TARGET_PRN,
        "methods": list(METHODS),
        "scenes": {key: list(value[0]) for key, value in SCENES.items()},
        "seeds": list(args.seeds),
        "truth_used_for_estimation": False,
        "truth_used_for_scoring_only": True,
        "observation_level": "window_average_401_tap",
        "notes": [
            "External baseline schedules were derived from one retained 401-point window average per second.",
            "This is a common-observable PVT result, but not a per-epoch baseline-estimator result.",
        ],
        "per_seed": results,
        "aggregate": aggregate,
    }
    output = args.output_root / "EXTERNAL_BASELINE_PVT_SUMMARY.json"
    output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"summary": str(output), "seeds": len(results)}, indent=2))


if __name__ == "__main__":
    main()
