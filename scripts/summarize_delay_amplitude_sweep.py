#!/usr/bin/env python3
"""Summarize delay-estimation accuracy over the delay/amplitude sweep."""

from __future__ import annotations

import csv
import math
import sys
from collections import defaultdict
from pathlib import Path


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit(f"usage: {sys.argv[0]} PER_CASE.csv SUMMARY.csv")
    source = Path(sys.argv[1])
    destination = Path(sys.argv[2])
    groups: dict[float, list[dict[str, str]]] = defaultdict(list)
    with source.open(newline="") as stream:
        for row in csv.DictReader(stream):
            groups[float(row["delta_tau_true"])].append(row)

    rows = []
    for delta in sorted(groups):
        group = groups[delta]
        valid = [row for row in group if row["scan_status"] == "OK"]
        errors = [float(row["delta_tau_hat"]) - delta for row in valid]
        estimates = [float(row["delta_tau_hat"]) for row in valid]
        mean = sum(estimates) / len(estimates) if estimates else float("nan")
        delay_rmse = math.sqrt(sum(error * error for error in errors) / len(errors)) if errors else float("nan")
        rows.append(
            {
                "delta_tau_true": f"{delta:.2f}",
                "n_valid": str(len(valid)),
                "n_failed": str(len(group) - len(valid)),
                "delta_tau_hat_mean": f"{mean:.5f}" if estimates else "",
                "delta_tau_hat_std": f"{math.sqrt(sum((estimate - mean) ** 2 for estimate in estimates) / len(estimates)):.5f}" if estimates else "",
                "delay_bias": f"{sum(errors) / len(errors):+.5f}" if errors else "",
                "delay_rmse": f"{delay_rmse:.5f}" if errors else "",
            }
        )

    fieldnames = [
        "delta_tau_true",
        "n_valid",
        "n_failed",
        "delta_tau_hat_mean",
        "delta_tau_hat_std",
        "delay_bias",
        "delay_rmse",
    ]
    with destination.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print("delta_tau_true  valid  failed  hat_mean  hat_std  bias      RMSE")
    for row in rows:
        print(
            f"{row['delta_tau_true']:>14}  {row['n_valid']:>5}  {row['n_failed']:>6}  "
            f"{row['delta_tau_hat_mean']:>8}  {row['delta_tau_hat_std']:>7}  "
            f"{row['delay_bias']:>8}  {row['delay_rmse']:>7}"
        )
    print(f"summary={destination}")


if __name__ == "__main__":
    main()
