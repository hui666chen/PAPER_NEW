#!/usr/bin/env python3
"""Run the frozen PRN28 observation-only gate on five new seeds."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np

ROOT = Path(__file__).resolve().parent
PAPER_ROOT = ROOT.parent
PRN30_ROOT = PAPER_ROOT / "gps_l1_prn30_shared_ratio_v2"
sys.path.insert(0, str(PAPER_ROOT))

TARGET_PRN = 28
SIMULATOR = Path("<LOCAL_USER>/chh/GPU/High/signalsim-new-multipath/build-if/IFdataGen")
RECEIVER = Path("<LOCAL_USER>/chh/GPU/High/gnss-sdr-cuda-scan/build-cuda-scan/src/main/gnss-sdr")
SCANNER = Path("<LOCAL_USER>/chh/GPU/High/gnss-sdr-cuda-scan/build-cuda-scan/src/algorithms/tracking/libs/frozen_reference_open_loop_scan")
CAL_CONFIG = ROOT / "configs/signalsim_prn28_isolated_145s.json"
VALID_CONFIG = ROOT / "configs/signalsim_prn28_validation_145s.json"
STATE_TEMPLATE = PAPER_ROOT / "gps_l1_prn30_formal_pvt/final_sign_persistence/state/run.conf"
BANKS = ROOT / "response_banks_25m"
CAL_AUDIT = ROOT / "CALIBRATION_AUDIT_25M.json"
THRESHOLD_FILE = PRN30_ROOT / "null_calibration/threshold_freeze.json"
RULE_FILE = PRN30_ROOT / "PERSISTENCE_RULE_FREEZE.json"
PROTOCOL = ROOT / "PRN28_MID_ELEVATION_PROTOCOL.md"
TEMP = Path("<TMP_DIR>/prn28_mid_elevation_batches")
CAL_IQ = TEMP / "prn28_isolated_current_binary_145s_iq16.bin"
CONFIRMATION_ROOT = ROOT / "confirmation_observation"
CONFIRMATION_SEEDS = (2026101120, 2026101121, 2026101122, 2026101123, 2026101124)
EXPECTED_IQ_BYTES = 14_499_900_000
SCENES = {
    "clean": (tuple(range(61, 77)), 68),
    "amp_0p05": (tuple(range(100, 116)), 108),
    "amp_0p10": (tuple(range(121, 137)), 129),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def command(command: list[str], log: Path, cwd: Path = ROOT) -> None:
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("w", encoding="utf-8", errors="replace") as stream:
        subprocess.run(command, cwd=cwd, stdout=stream,
                       stderr=subprocess.STDOUT, check=True)


def replace_or_append(lines: list[str], key: str, value: str) -> None:
    prefix = key + "="
    for index, line in enumerate(lines):
        if line.startswith(prefix):
            lines[index] = prefix + value
            return
    lines.append(prefix + value)


def generate(config: Path, output: Path, seed: int, log: Path) -> None:
    if output.exists():
        raise RuntimeError(f"refusing to overwrite {output}")
    command([str(SIMULATOR), "--config", str(config), "--output", str(output),
             "--noise-seed", str(seed), "--multi-thread"], log)
    if output.stat().st_size != EXPECTED_IQ_BYTES:
        raise RuntimeError(f"unexpected IQ size {output.stat().st_size}")
    if "Clipped samples: 0" not in log.read_text(errors="replace"):
        raise RuntimeError(f"{config}: simulator clipping check failed")


def prepare_state(run_dir: Path, iq: Path) -> Path:
    state_dir = run_dir / "state"
    state_dir.mkdir(parents=True, exist_ok=True)
    lines = STATE_TEMPLATE.read_text(encoding="utf-8").splitlines()
    replace_or_append(lines, "SignalSource.filename", str(iq))
    replace_or_append(lines, "Channel0.satellite", str(TARGET_PRN))
    replace_or_append(lines, "Tracking_1C.corr_scan_prn", str(TARGET_PRN))
    replace_or_append(lines, "Tracking_1C.state_dump_filename",
                      str(state_dir / "tracking_state_prn28.csv"))
    replace_or_append(lines, "Observables.dump", "false")
    replace_or_append(lines, "PVT.dump", "false")
    path = run_dir / "configs/state.conf"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def extract_scene(run_dir: Path, scene: str, windows: tuple[int, ...],
                  validation_iq: Path, state: Path) -> None:
    destination = run_dir / "open_loop" / scene
    destination.mkdir(parents=True, exist_ok=True)
    output = destination / "prn28_epoch_responses.csv"
    command([
        str(SCANNER), "--prn", str(TARGET_PRN), "--state", str(state),
        "--primary", str(CAL_IQ), "--secondary", str(CAL_IQ),
        "--mixed", str(validation_iq), "--output-dir", str(destination),
        "--epoch-output", str(output), "--epochs-per-window", "100",
        "--max-windows", "16", "--grid-ms", "1000",
        "--skip-ms", str(windows[0] * 1000), "--fs-hz", "25000000",
        "--if-hz", "-4580000",
    ], destination / "scan.log")
    expected = 16 * 100 * 401 + 1
    rows = sum(1 for _ in output.open(encoding="utf-8"))
    if rows != expected:
        raise RuntimeError(f"{scene}: {rows} rows, expected {expected}")


def run_diagnostic(run_dir: Path) -> Path:
    output = run_dir / "shared_ratio"
    command([sys.executable, str(ROOT / "shared_ratio_existence_prn28.py"),
             "--run-dir", str(run_dir), "--banks-dir", str(BANKS),
             "--output", str(output)], run_dir / "shared_ratio.log")
    return output


def apply_rule(records: list[dict], k: int, m: int) -> list[dict]:
    output = [dict(row) for row in records]
    for scene in SCENES:
        rows = sorted((row for row in output if row["scene"] == scene),
                      key=lambda row: int(row["window"]))
        evidence = [bool(row["single_window_exists"]) for row in rows]
        for index, row in enumerate(rows):
            active = index >= m - 1 and sum(evidence[index - m + 1:index + 1]) >= k
            row["gate_active"] = bool(active)
            row["rule_k"] = int(k)
            row["rule_m"] = int(m)
    return output


def summarize(records: list[dict]) -> dict:
    result = {}
    for scene in SCENES:
        rows = [row for row in records if row["scene"] == scene]
        result[scene] = {
            "windows": len(rows),
            "single_window_exists": int(sum(row["single_window_exists"] for row in rows)),
            "gate_active": int(sum(row["gate_active"] for row in rows)),
            "gate_active_rate": float(sum(row["gate_active"] for row in rows) / len(rows)),
            "cv_no_boundary": int(sum(row["cv_no_boundary"] for row in rows)),
            "cv_gain_min_quantiles": [float(value) for value in np.quantile(
                [row["cv_gain_min"] for row in rows], [0.1, 0.5, 0.9]
            )],
        }
    return result


def process_seed(seed: int, threshold: float, k: int, m: int) -> dict:
    run_dir = CONFIRMATION_ROOT / f"seed_{seed}"
    if run_dir.exists():
        raise RuntimeError(f"refusing to overwrite {run_dir}")
    run_dir.mkdir(parents=True)
    validation_iq = TEMP / f"validation_seed_{seed}_145s_iq16.bin"
    try:
        print(f"seed={seed} stage=generate", flush=True)
        generate(VALID_CONFIG, validation_iq, seed, run_dir / "generation.log")
        print(f"seed={seed} stage=state", flush=True)
        state_config = prepare_state(run_dir, validation_iq)
        command([str(RECEIVER), f"--config_file={state_config}",
                 "--logtostderr=1", "--colorlogtostderr=1"],
                run_dir / "state_console.log", state_config.parent)
        state = run_dir / "state/tracking_state_prn28.csv"
        if not state.is_file() or state.stat().st_size == 0:
            raise RuntimeError(f"missing tracking state {state}")
        print(f"seed={seed} stage=responses", flush=True)
        for scene, (windows, _) in SCENES.items():
            extract_scene(run_dir, scene, windows, validation_iq, state)
        print(f"seed={seed} stage=shared_ratio", flush=True)
        output = run_diagnostic(run_dir)
        records = json.loads((output / "records.json").read_text())
        for row in records:
            row["single_window_exists"] = bool(
                row["cv_gain_min"] > threshold and row["cv_no_boundary"]
            )
        records = apply_rule(records, k, m)
        (output / "gated_records.json").write_text(
            json.dumps(records, indent=2) + "\n", encoding="utf-8"
        )
        summary = {
            "seed": int(seed), "truth_accessed": False,
            "threshold": float(threshold), "rule": {"k": int(k), "m": int(m)},
            "state_sha256": sha256(state), "scenes": summarize(records),
        }
        (run_dir / "persistence_summary.json").write_text(
            json.dumps(summary, indent=2) + "\n", encoding="utf-8"
        )
        for scene in SCENES:
            (run_dir / "open_loop" / scene / "prn28_epoch_responses.csv").unlink()
        return summary
    finally:
        if validation_iq.exists() and validation_iq.parent == TEMP:
            validation_iq.unlink()


def main() -> None:
    for path in (SIMULATOR, RECEIVER, SCANNER, CAL_CONFIG, VALID_CONFIG,
                 STATE_TEMPLATE, BANKS, CAL_AUDIT, THRESHOLD_FILE, RULE_FILE,
                 PROTOCOL):
        if not path.exists():
            raise FileNotFoundError(path)
    if CONFIRMATION_ROOT.exists():
        raise RuntimeError(f"refusing to overwrite {CONFIRMATION_ROOT}")
    audit = json.loads(CAL_AUDIT.read_text())
    if audit["target_prn"] != TARGET_PRN:
        raise RuntimeError("PRN28 calibration target mismatch")
    if audit["receiver_sha256"] != sha256(RECEIVER):
        raise RuntimeError("current receiver binary differs from calibration audit")
    rule = json.loads(RULE_FILE.read_text())["selection"]["selected"]
    k, m = int(rule["k"]), int(rule["m"])
    if (k, m) != (3, 5):
        raise RuntimeError(f"unexpected frozen rule {(k, m)}")
    threshold = float(json.loads(THRESHOLD_FILE.read_text())["T_null"])
    TEMP.mkdir(parents=True, exist_ok=True)
    if any(TEMP.iterdir()):
        raise RuntimeError(f"temporary directory is not empty: {TEMP}")
    CONFIRMATION_ROOT.mkdir(parents=True)
    summaries = []
    try:
        print("stage=generate_calibration_primary", flush=True)
        generate(CAL_CONFIG, CAL_IQ, int(audit["calibration_seed"]),
                 TEMP / "confirmation_calibration_primary_generation.log")
        if sha256(CAL_IQ) != audit["calibration_iq_sha256"]:
            raise RuntimeError("regenerated calibration IQ hash differs from audit")
        for seed in CONFIRMATION_SEEDS:
            if shutil.disk_usage("/tmp").free < 18_000_000_000:
                raise RuntimeError("less than 18 GB available for next seed")
            summaries.append(process_seed(seed, threshold, k, m))
        payload = {
            "status": "prn28_observation_gate_confirmation_complete",
            "truth_accessed": False, "target_prn": TARGET_PRN,
            "calibration_audit_sha256": sha256(CAL_AUDIT),
            "threshold_file_sha256": sha256(THRESHOLD_FILE),
            "rule_file_sha256": sha256(RULE_FILE),
            "protocol_sha256": sha256(PROTOCOL), "threshold": threshold,
            "rule": {"k": k, "m": m}, "seeds": list(CONFIRMATION_SEEDS),
            "seed_summaries": summaries,
        }
        (CONFIRMATION_ROOT / "CONFIRMATION_SUMMARY.json").write_text(
            json.dumps(payload, indent=2) + "\n", encoding="utf-8"
        )
        print(json.dumps({"rule": payload["rule"],
                          "summary": str(CONFIRMATION_ROOT / "CONFIRMATION_SUMMARY.json")},
                         indent=2), flush=True)
    finally:
        if CAL_IQ.exists() and CAL_IQ.parent == TEMP:
            CAL_IQ.unlink()
        if TEMP.exists() and not any(TEMP.iterdir()):
            TEMP.rmdir()


if __name__ == "__main__":
    main()
