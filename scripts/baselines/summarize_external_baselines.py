#!/usr/bin/env python3
"""Build compact schedules and summaries for the external baselines.

The input CSVs contain only receiver observations and estimator outputs.  This
script does not read simulator truth.  It converts each method's epoch/window
delays into the frozen one-second correction schedule used by the common-PVT
adapter and records the boundary/validity audit separately.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np


CHIP_METRES = 299_792_458.0 / 1.023e6
METHOD_COLUMNS = {
    "eml": ("eml_tau_chip", "eml_boundary"),
    "double_delta": ("dd_tau_chip", "dd_boundary"),
    "medll_ideal": ("medll_ideal_tau_chip", "medll_ideal_boundary"),
    "medll_chain": ("medll_chain_tau_chip", "medll_chain_boundary"),
}
SCENES = ("clean", "amp_0p05", "amp_0p10")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_schedule(
    path: Path, rows: list[dict[str, str]], method: str, prn: int, sample_rate_hz: int,
    min_valid_fraction: float, offset_chip: float,
) -> dict[str, object]:
    tau_key, boundary_key = METHOD_COLUMNS[method]
    groups: dict[int, list[dict[str, str]]] = {}
    for row in rows:
        groups.setdefault(int(row["window"]), []).append(row)

    output: list[dict[str, object]] = []
    accepted_windows = 0
    boundary_count = 0
    tau_values: list[float] = []
    for window in sorted(groups):
        group = groups[window]
        boundary_count += sum(
            row[boundary_key].strip().lower() in {"1", "true", "yes"}
            for row in group
        )
        valid = []
        for row in group:
            try:
                tau = float(row[tau_key])
            except (KeyError, ValueError):
                continue
            if row[boundary_key].strip().lower() in {"1", "true", "yes"}:
                continue
            if np.isfinite(tau):
                valid.append(tau)
        accepted = len(valid) >= min_valid_fraction * len(group)
        tau = float(np.median(valid)) if accepted else 0.0
        relative_tau = tau - offset_chip if accepted else 0.0
        if accepted:
            accepted_windows += 1
            tau_values.append(tau)
        output.append({
            "start_sample": window * sample_rate_hz,
            "end_sample": (window + 1) * sample_rate_hz - 1,
            "prn": prn,
            "correction_m": f"{-relative_tau * CHIP_METRES:.12f}",
            "valid_epoch_count": len(valid),
            "epoch_count": len(group),
            "accepted": accepted,
        })

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(output[0]))
        writer.writeheader()
        writer.writerows(output)
    return {
        "windows": len(output),
        "accepted_windows": accepted_windows,
        "rejected_windows": len(output) - accepted_windows,
        "boundary_rate": boundary_count / max(len(rows), 1),
        "tau_median_chip": float(np.median(tau_values)) if tau_values else None,
        "relative_tau_median_chip": (
            float(np.median(np.asarray(tau_values) - offset_chip))
            if tau_values else None
        ),
        "offset_chip": offset_chip,
        "correction_abs_median_m": (
            float(np.median(np.abs((np.asarray(tau_values) - offset_chip) * CHIP_METRES)))
            if tau_values else None
        ),
        "rows": len(rows),
        "rows_per_window_min": min(len(group) for group in groups.values()),
        "rows_per_window_max": max(len(group) for group in groups.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--prn", type=int, default=9)
    parser.add_argument("--sample-rate-hz", type=int, default=25_000_000)
    parser.add_argument("--min-valid-fraction", type=float, default=0.5)
    args = parser.parse_args()
    if not 0.0 < args.min_valid_fraction <= 1.0:
        raise ValueError("min-valid-fraction must be in (0, 1]")

    schedules = args.output_root / "schedules"
    records: list[dict[str, object]] = []
    inputs: dict[tuple[int, str], Path] = {}
    for input_path in sorted(args.input_root.glob("seed_*_*.csv")):
        stem = input_path.stem
        parts = stem.split("_")
        if len(parts) < 3 or parts[0] != "seed":
            continue
        seed = int(parts[1])
        scene = "_".join(parts[2:])
        if scene not in SCENES:
            continue
        inputs[(seed, scene)] = input_path

    for (seed, scene), input_path in sorted(inputs.items()):
        with input_path.open(newline="", encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
        clean_path = inputs.get((seed, "clean"))
        if clean_path is None:
            raise ValueError(f"seed {seed}: clean reference is required")
        with clean_path.open(newline="", encoding="utf-8") as stream:
            clean_rows = list(csv.DictReader(stream))
        for method in METHOD_COLUMNS:
            tau_key, boundary_key = METHOD_COLUMNS[method]
            clean_valid = [
                float(row[tau_key]) for row in clean_rows
                if row[boundary_key].strip().lower() not in {"1", "true", "yes"}
                and np.isfinite(float(row[tau_key]))
            ]
            if not clean_valid:
                raise ValueError(f"seed {seed}/{method}: no valid clean reference")
            offset_chip = float(np.median(clean_valid))
            schedule_path = schedules / f"seed_{seed}_{scene}_{method}.csv"
            audit = write_schedule(
                schedule_path, rows, method, args.prn, args.sample_rate_hz,
                args.min_valid_fraction, offset_chip,
            )
            records.append({
                "seed": seed,
                "scene": scene,
                "method": method,
                "input_csv": str(input_path),
                "input_sha256": sha256(input_path),
                "schedule_csv": str(schedule_path),
                "schedule_sha256": sha256(schedule_path),
                "observation_level": "window_average_401_tap" if audit["rows_per_window_min"] == 1 else "epoch_401_tap",
                "zero_reference": str(clean_path),
                "zero_reference_sha256": sha256(clean_path),
                **audit,
            })

    args.output_root.mkdir(parents=True, exist_ok=True)
    manifest = {
        "status": "external_baseline_observation_summaries_complete",
        "truth_accessed": False,
        "pvt_scored": False,
        "prn": args.prn,
        "chip_metres": CHIP_METRES,
        "min_valid_fraction": args.min_valid_fraction,
        "methods": list(METHOD_COLUMNS),
        "records": records,
    }
    output = args.output_root / "EXTERNAL_BASELINE_SUMMARY.json"
    output.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"summary": str(output), "records": len(records)}, indent=2))


if __name__ == "__main__":
    main()
