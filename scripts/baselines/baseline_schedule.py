#!/usr/bin/env python3
"""Convert observation-level baseline delays into a PVT correction schedule.

The schedule rule is frozen for the baseline comparison: use the median valid
epoch estimate per one-second window; if fewer than half the epochs are valid,
write zero correction and report the window as rejected.  This is a
deterministic validity rule, not the proposed 3-of-5 multipath gate.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

CHIP_METRES = 299_792_458.0 / 1.023e6


METHOD_COLUMNS = {
    "eml": ("eml_tau_chip", "eml_boundary"),
    "double_delta": ("dd_tau_chip", "dd_boundary"),
    "medll_ideal": ("medll_ideal_tau_chip", "medll_ideal_boundary"),
    "medll_chain": ("medll_chain_tau_chip", "medll_chain_boundary"),
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-csv", type=Path, required=True)
    parser.add_argument("--output-csv", type=Path, required=True)
    parser.add_argument("--method", choices=sorted(METHOD_COLUMNS), required=True)
    parser.add_argument("--prn", type=int, required=True)
    parser.add_argument("--sample-rate-hz", type=int, default=25_000_000)
    parser.add_argument("--min-valid-fraction", type=float, default=0.5)
    args = parser.parse_args()
    if not 0.0 < args.min_valid_fraction <= 1.0:
        raise ValueError("min-valid-fraction must be in (0, 1]")

    tau_key, boundary_key = METHOD_COLUMNS[args.method]
    with args.baseline_csv.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    groups: dict[int, list[dict[str, str]]] = {}
    for row in rows:
        groups.setdefault(int(row["window"]), []).append(row)
    output = []
    rejected = 0
    for window in sorted(groups):
        group = groups[window]
        valid = [
            float(row[tau_key]) for row in group
            if row[boundary_key].lower() in {"0", "false", "no"}
            and np.isfinite(float(row[tau_key]))
        ]
        accepted = len(valid) >= args.min_valid_fraction * len(group)
        tau = float(np.median(valid)) if accepted else 0.0
        if not accepted:
            rejected += 1
        start = window * args.sample_rate_hz
        output.append({
            "start_sample": start,
            "end_sample": (window + 1) * args.sample_rate_hz - 1,
            "prn": args.prn,
            "correction_m": -tau * CHIP_METRES,
            "valid_epoch_count": len(valid),
            "epoch_count": len(group),
            "accepted": accepted,
        })
    if not output:
        raise ValueError("baseline CSV contains no rows")
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_csv.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(output[0]))
        writer.writeheader()
        writer.writerows(output)
    print({"method": args.method, "windows": len(output), "rejected_windows": rejected})


if __name__ == "__main__":
    main()
