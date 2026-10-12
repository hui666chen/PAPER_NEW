#!/usr/bin/env python3
"""Benchmark observation-level baseline cost and record receiver RT factor.

The script intentionally separates algorithmic benchmark time from the
end-to-end CUDA receiver time.  The latter is entered from an unfiltered
receiver log or a timed run, never inferred from a Python prototype.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

import numpy as np

from correlator_baselines import (
    BaselineConfig,
    EpochScan,
    IdealTemplate,
    eml_estimate,
    double_delta_estimate,
    medll_estimate,
)


def synthetic_scan(taps: int = 401) -> EpochScan:
    x = np.linspace(-2.0, 2.0, taps)
    primary = np.maximum(1.0 - np.abs(x - 0.07), 0.0)
    reflected = 0.12 * np.exp(1j * 0.7) * np.maximum(1.0 - np.abs(x + 0.43), 0.0)
    noise = np.random.default_rng(20261012).normal(size=taps) * 1.0e-3
    return EpochScan(0, 0, 0, x, primary + reflected + noise, {})


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--repeat", type=int, default=20)
    parser.add_argument("--receiver-wall-seconds", type=float)
    parser.add_argument("--signal-seconds", type=float, default=145.0)
    args = parser.parse_args()
    if args.repeat < 1:
        raise ValueError("repeat must be positive")

    cfg = BaselineConfig()
    scan = synthetic_scan()
    template = IdealTemplate(scan.delay_chip, 0.07)
    timings = {}
    for name, function in (
        ("eml_0p10", lambda: eml_estimate(scan, cfg)),
        ("double_delta_0p10_0p50", lambda: double_delta_estimate(scan, cfg)),
        ("medll_ideal", lambda: medll_estimate(scan, template, cfg)),
    ):
        started = time.perf_counter()
        for _ in range(args.repeat):
            function()
        elapsed = time.perf_counter() - started
        timings[name] = {
            "repeat": args.repeat,
            "total_seconds": elapsed,
            "seconds_per_epoch": elapsed / args.repeat,
            "implementation": "Python reference, not the CUDA receiver",
        }
    receiver = None
    if args.receiver_wall_seconds is not None:
        receiver = {
            "wall_seconds": args.receiver_wall_seconds,
            "signal_seconds": args.signal_seconds,
            "real_time_factor": args.signal_seconds / args.receiver_wall_seconds,
            "definition": "processed signal duration divided by wall-clock receiver runtime",
        }
    payload = {
        "truth_accessed": False,
        "algorithmic_reference_timings": timings,
        "end_to_end_cuda_receiver": receiver,
        "warning": "Do not report the Python reference timings as GPU timing.",
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
