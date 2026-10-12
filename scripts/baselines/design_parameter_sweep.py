#!/usr/bin/env python3
"""Create the frozen multipath-baseline parameter design."""

from __future__ import annotations

import argparse
import csv
from itertools import product
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--delays", default="0.05,0.10,0.20,0.40,0.60,0.80,1.00"
    )
    parser.add_argument("--ratios", default="0.10,0.20,0.34,0.50,0.70")
    parser.add_argument("--phases", default="0,90,180,270")
    parser.add_argument("--cn0", default="35,40,45,50")
    args = parser.parse_args()

    values = {
        "delay_chip": [float(value) for value in args.delays.split(",") if value],
        "amplitude_ratio": [float(value) for value in args.ratios.split(",") if value],
        "phase_deg": [float(value) for value in args.phases.split(",") if value],
        "cn0_db_hz": [float(value) for value in args.cn0.split(",") if value],
    }
    rows = [dict(zip(values, item)) for item in product(*values.values())]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(values))
        writer.writeheader()
        writer.writerows(rows)
    print(f"rows={len(rows)} output={args.output}")


if __name__ == "__main__":
    main()
