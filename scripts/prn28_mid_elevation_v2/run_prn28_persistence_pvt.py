#!/usr/bin/env python3
"""Run common-observable three-way PVT for the PRN28 middle-elevation batch."""

from __future__ import annotations

import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
PAPER_ROOT = ROOT.parent
PRN30_ROOT = PAPER_ROOT / "gps_l1_prn30_shared_ratio_v2"
FORMAL_ROOT = PAPER_ROOT / "gps_l1_prn30_formal_pvt"
sys.path.insert(0, str(PRN30_ROOT))
sys.path.insert(0, str(PAPER_ROOT))

import run_prn30_persistence_pvt as persistence  # noqa: E402
import run_prn30_v2_batch as batch  # noqa: E402
import run_prn30_v2_pvt as base  # noqa: E402


PVT_ROOT = ROOT / "persistence_pvt_confirmation"
CONFIRM_ROOT = ROOT / "confirmation_observation"
RULE_FILE = PRN30_ROOT / "PERSISTENCE_RULE_FREEZE.json"
THRESHOLD_FILE = PRN30_ROOT / "null_calibration/threshold_freeze.json"
PROTOCOL = ROOT / "PRN28_MID_ELEVATION_PROTOCOL.md"
CAL_AUDIT = ROOT / "CALIBRATION_AUDIT_25M.json"
VALID_CONFIG = ROOT / "configs/signalsim_prn28_validation_145s.json"
RECEIVER_TEMPLATE = ROOT / "configs/receiver_prn28_source.conf"
STATE_TEMPLATE = FORMAL_ROOT / "final_sign_persistence/state/run.conf"
CONFIRM_SEEDS = (2026101120, 2026101121, 2026101122, 2026101123, 2026101124)
METHODS = ("single", "always_shared_ratio", "gated_persistence_shared_ratio")
SCENES = {
    "clean": (tuple(range(61, 77)), 68),
    "amp_0p05": (tuple(range(100, 116)), 108),
    "amp_0p10": (tuple(range(121, 137)), 129),
}
EXPECTED_SATELLITES = (
    "G05", "G07", "G08", "G09", "G13", "G14", "G17", "G20", "G21", "G28", "G30"
)
TEMP = Path("<TMP_DIR>/prn28_mid_elevation_pvt")


def configure_imports() -> None:
    base.TARGET_PRN = 28
    base.TARGET_SATELLITE = "G28"
    base.METHODS = METHODS
    base.SCENES = SCENES
    base.CONFIRMATION_SEEDS = CONFIRM_SEEDS
    base.EXPECTED_SATELLITES = EXPECTED_SATELLITES
    base.EXPECTED_SATELLITE_COUNT = len(EXPECTED_SATELLITES)
    base.DYNAMIC_EXPECTED_SATELLITES = False
    base.PVT_ROOT = PVT_ROOT
    base.CONFIRM_ROOT = CONFIRM_ROOT
    base.PROTOCOL = PROTOCOL
    base.THRESHOLD_FILE = THRESHOLD_FILE
    base.CAL_AUDIT = CAL_AUDIT
    base.RECEIVER_TEMPLATE = RECEIVER_TEMPLATE
    base.STATE_TEMPLATE = STATE_TEMPLATE
    base.TEMP = TEMP

    batch.TARGET_PRN = 28
    batch.SIMULATOR = Path("<LOCAL_USER>/chh/GPU/High/signalsim-new-multipath/build-if/IFdataGen")
    batch.RECEIVER = Path("<LOCAL_USER>/chh/GPU/High/gnss-sdr-cuda-scan/build-cuda-scan/src/main/gnss-sdr")
    batch.VALID_CONFIG = VALID_CONFIG
    batch.STATE_TEMPLATE = STATE_TEMPLATE
    batch.EXPECTED_IQ_BYTES = 14_499_900_000
    batch.SCENES = SCENES
    batch.TEMP = TEMP
    batch.CAL_AUDIT = CAL_AUDIT

    persistence.PVT_ROOT = PVT_ROOT
    persistence.CONFIRM_ROOT = CONFIRM_ROOT
    persistence.RULE_FILE = RULE_FILE
    persistence.THRESHOLD_FILE = THRESHOLD_FILE
    persistence.PROTOCOL = PROTOCOL
    persistence.CONFIRM_SEEDS = CONFIRM_SEEDS
    persistence.METHODS = METHODS
    persistence.TEMP = TEMP
    persistence.base = base
    persistence.batch = batch


def main() -> None:
    configure_imports()
    persistence.main()
    summary_path = PVT_ROOT / "PVT_CONFIRMATION_SUMMARY.json"
    payload = json.loads(summary_path.read_text())
    payload["status"] = "prn28_mid_elevation_persistence_pvt_confirmation_complete"
    payload["reference_pipeline"] = "run_prn30_persistence_pvt with PRN28 paths and 11-satellite geometry"
    summary_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"PRN28_COMPLETE={PVT_ROOT}", flush=True)


if __name__ == "__main__":
    main()
