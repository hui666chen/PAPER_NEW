#!/usr/bin/env python3
"""Observation-only shared-ratio diagnostic for PRN28 validation."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parent
PAPER_ROOT = ROOT.parent
sys.path.insert(0, str(PAPER_ROOT))
sys.path.insert(0, str(PAPER_ROOT / "gps_l1_prn09_comparable_v3"))

import shared_ratio_existence_prn09_comparable as base  # noqa: E402
from receiver_chain_response_control import ResponseSpec, load_response_model  # noqa: E402
from monte_carlo_resolution_baseline import build_dictionary  # noqa: E402

TARGET_PRN = 28
SPEC = ResponseSpec(prn=TARGET_PRN)
SCENES = {
    "clean": (tuple(range(61, 77)), 68),
    "amp_0p05": (tuple(range(100, 116)), 108),
    "amp_0p10": (tuple(range(121, 137)), 129),
}

base.TARGET_PRN = TARGET_PRN
base.SPEC = SPEC
base.SCENES = SCENES
base.base.TARGET_PRN = TARGET_PRN
base.base.SPEC = SPEC
base.base.SCENES = SCENES


def run(run_dir: Path, banks_dir: Path, output: Path) -> dict:
    if output.exists():
        raise RuntimeError(f"refusing to overwrite {output}")
    output.mkdir(parents=True)
    rho_values = base.base.make_rho_grid()
    records = []
    for scene, (windows, bank_window) in SCENES.items():
        source = run_dir / "open_loop" / scene / "prn28_epoch_responses.csv"
        observations = base.base.load_windows(source, windows)
        model = load_response_model(
            banks_dir / f"bank_{bank_window}.npz", expected_spec=SPEC
        )
        delay_grid, dictionary, _, _, _ = build_dictionary(model, SPEC)
        for window in windows:
            values, weights, magnitude = observations[window]
            record = base.base.analyze_window(
                values, weights, dictionary, delay_grid, rho_values
            )
            record.update({
                "scene": scene,
                "window": int(window),
                "bank_window": int(bank_window),
                "prompt_magnitude_median": float(np.median(magnitude)),
            })
            records.append(record)
            print(json.dumps(record), flush=True)

    records.sort(key=lambda row: (row["scene"], row["window"]))
    (output / "records.json").write_text(
        json.dumps(records, indent=2) + "\n", encoding="utf-8"
    )
    with (output / "records.csv").open("w", newline="", encoding="utf-8") as stream:
        fields = list(records[0])
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for row in records:
            flattened = dict(row)
            flattened["A_to_B"] = json.dumps(row["A_to_B"], sort_keys=True)
            flattened["B_to_A"] = json.dumps(row["B_to_A"], sort_keys=True)
            writer.writerow(flattened)

    summary = {}
    for scene, (windows, _) in SCENES.items():
        rows = [row for row in records if row["scene"] == scene]
        summary[scene] = {
            "windows": len(rows),
            "cv_both_positive_count": int(sum(row["cv_both_positive"] for row in rows)),
            "cv_no_boundary_count": int(sum(row["cv_no_boundary"] for row in rows)),
            "cv_gain_min_quantiles": [
                float(v) for v in np.quantile(
                    [row["cv_gain_min"] for row in rows], [0.1, 0.5, 0.9]
                )
            ],
            "rho_magnitude_mean_quantiles": [
                float(v) for v in np.quantile(
                    [row["rho_magnitude_mean"] for row in rows], [0.1, 0.5, 0.9]
                )
            ],
        }
    payload = {
        "target_prn": TARGET_PRN,
        "truth_accessed": False,
        "purpose": "observation-only shared-ratio diagnostic for frozen PVT gate",
        "scene_subset": list(SCENES),
        "split": "even epochs train/odd epochs heldout and swapped",
        "model": "per-epoch common gain and nuisance coefficients; shared complex reflection/direct ratio",
        "rho_grid": [
            {"real": float(z.real), "imag": float(z.imag), "magnitude": float(abs(z))}
            for z in rho_values
        ],
        "summary": summary,
    }
    (output / "summary.json").write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8"
    )
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--banks-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    payload = run(args.run_dir.resolve(), args.banks_dir.resolve(), args.output.resolve())
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
