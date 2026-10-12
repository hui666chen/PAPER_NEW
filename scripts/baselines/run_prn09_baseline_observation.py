#!/usr/bin/env python3
"""Regenerate a PRN9 observation batch and retain mixed 401-tap CSVs.

This wrapper reuses the already frozen PRN9 receiver/simulator/calibration
pipeline.  It is intentionally separate from the published confirmation
batch: retaining large response CSVs is opt-in and the output root must be a
new empty directory.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1] / "prn09_comparable_v3"
sys.path.insert(0, str(ROOT))
import run_prn09_comparable_batch as batch  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--temp-root", type=Path, required=True)
    parser.add_argument(
        "--experiment-root", type=Path,
        default=Path("/home/zoxoai/Documents/ChatGPT/paper/gps_l1_prn09_comparable_v3"),
        help="local PRN9 experiment directory containing configs, banks, and audit files",
    )
    parser.add_argument(
        "--formal-pvt-root", type=Path,
        default=Path("/home/zoxoai/Documents/ChatGPT/paper/gps_l1_prn30_formal_pvt"),
    )
    parser.add_argument(
        "--shared-ratio-root", type=Path,
        default=Path("/home/zoxoai/Documents/ChatGPT/paper/gps_l1_prn30_shared_ratio_v2"),
    )
    parser.add_argument(
        "--high-root", type=Path,
        default=Path("/home/zoxoai/chh/GPU/High"),
        help="local GPU receiver/simulator build root",
    )
    parser.add_argument(
        "--seeds", nargs="+", type=int,
        default=[2026101110, 2026101111, 2026101112, 2026101113, 2026101114],
    )
    args = parser.parse_args()
    if args.output_root.exists() and any(args.output_root.iterdir()):
        raise RuntimeError(f"refusing non-empty output root: {args.output_root}")
    if args.temp_root.exists() and any(args.temp_root.iterdir()):
        raise RuntimeError(f"refusing non-empty temp root: {args.temp_root}")

    batch.CONFIRMATION_ROOT = args.output_root
    batch.TEMP = args.temp_root
    batch.CAL_IQ = args.temp_root / "prn09_isolated_current_binary_145s_iq16.bin"
    batch.CONFIRMATION_SEEDS = tuple(args.seeds)
    batch.CAL_CONFIG = args.experiment_root / "configs/signalsim_prn09_isolated_145s.json"
    batch.VALID_CONFIG = args.experiment_root / "configs/signalsim_prn09_validation_145s.json"
    batch.BANKS = args.experiment_root / "response_banks"
    batch.CAL_AUDIT = args.experiment_root / "CALIBRATION_AUDIT.json"
    batch.PROTOCOL = args.experiment_root / "COMPARABILITY_PROTOCOL.md"
    batch.STATE_TEMPLATE = (
        args.formal_pvt_root / "final_sign_persistence/state/run.conf"
    )
    batch.THRESHOLD_FILE = (
        args.shared_ratio_root / "null_calibration/threshold_freeze.json"
    )
    batch.RULE_FILE = args.shared_ratio_root / "PERSISTENCE_RULE_FREEZE.json"
    batch.SIMULATOR = args.high_root / "signalsim-new-multipath/build-if/IFdataGen"
    batch.RECEIVER = (
        args.high_root / "gnss-sdr-cuda-scan/build-cuda-scan/src/main/gnss-sdr"
    )
    batch.SCANNER = (
        args.high_root / "gnss-sdr-cuda-scan/build-cuda-scan/src/algorithms/"
        "tracking/libs/frozen_reference_open_loop_scan"
    )
    batch.main(["--keep-epoch-responses"])
    print(f"BASELINE_OBSERVATIONS={args.output_root}")


if __name__ == "__main__":
    main()
