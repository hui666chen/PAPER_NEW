#!/usr/bin/env python3
"""Calibrate PRN28 response banks with the current receiver binary."""

from __future__ import annotations

import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
PAPER_ROOT = ROOT.parent
sys.path.insert(0, str(PAPER_ROOT))
sys.path.insert(0, str(PAPER_ROOT / "gps_l1_prn30_shared_ratio_v2"))

import run_prn30_v2_calibration as base  # noqa: E402
from receiver_chain_response_control import ResponseSpec  # noqa: E402


TARGET = 28
PROTOCOL = ROOT / "PRN28_MID_ELEVATION_PROTOCOL.md"
STATE_TEMPLATE = PAPER_ROOT / "gps_l1_prn30_formal_pvt/final_sign_persistence/state/run.conf"

base.ROOT = ROOT
base.TARGET_PRN = TARGET
base.SPEC = ResponseSpec(prn=TARGET)
base.CAL_CONFIG = ROOT / "configs/signalsim_prn28_isolated_145s.json"
base.STATE_TEMPLATE = STATE_TEMPLATE
base.TEMP = Path("<TMP_DIR>/prn28_mid_elevation_calibration")
base.CAL_IQ = base.TEMP / "prn28_isolated_current_binary_145s_iq16.bin"
base.OUT = ROOT / "calibration_25m"
base.BANKS = ROOT / "response_banks_25m"
base.PROTOCOL = PROTOCOL
base.MANIFEST = ROOT / "calibration_manifest.json"
base.CAL_SEED = 2026101118


def main() -> None:
    base.main()
    audit_path = ROOT / "CALIBRATION_AUDIT_25M.json"
    audit = json.loads(audit_path.read_text())
    audit["status"] = "current_binary_prn28_calibration_complete"
    audit["target_prn"] = TARGET
    audit["old_prn28_response_banks_reused"] = False
    audit["protocol_sha256"] = base.sha256(PROTOCOL)
    audit_path.write_text(json.dumps(audit, indent=2) + "\n")
    print(f"COMPLETE={audit_path}", flush=True)


if __name__ == "__main__":
    main()
