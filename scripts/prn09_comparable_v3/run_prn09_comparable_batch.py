#!/usr/bin/env python3
"""Run the frozen PRN9 comparable observation-only confirmation batch."""

from __future__ import annotations

import hashlib
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np

ROOT = Path(__file__).resolve().parent
PAPER_ROOT = ROOT.parent
sys.path.insert(0, str(PAPER_ROOT))

TARGET_PRN = 9
SIMULATOR = Path(
    "<LOCAL_USER>/chh/GPU/High/signalsim-new-multipath/build-if/IFdataGen"
)
RECEIVER = Path(
    "<LOCAL_USER>/chh/GPU/High/gnss-sdr-cuda-scan/"
    "build-cuda-scan/src/main/gnss-sdr"
)
SCANNER = Path(
    "<LOCAL_USER>/chh/GPU/High/gnss-sdr-cuda-scan/"
    "build-cuda-scan/src/algorithms/tracking/libs/"
    "frozen_reference_open_loop_scan"
)
CAL_CONFIG = ROOT / "configs/signalsim_prn09_isolated_145s.json"
VALID_CONFIG = ROOT / "configs/signalsim_prn09_validation_145s.json"
STATE_TEMPLATE = (
    PAPER_ROOT / "gps_l1_prn30_formal_pvt/final_sign_persistence/state/run.conf"
)
BANKS = ROOT / "response_banks"
CAL_AUDIT = ROOT / "CALIBRATION_AUDIT.json"
THRESHOLD_FILE = (
    PAPER_ROOT / "gps_l1_prn30_shared_ratio_v2/null_calibration/threshold_freeze.json"
)
RULE_FILE = PAPER_ROOT / "gps_l1_prn30_shared_ratio_v2/PERSISTENCE_RULE_FREEZE.json"
PROTOCOL = ROOT / "COMPARABILITY_PROTOCOL.md"
CONFIRMATION_ROOT = ROOT / "confirmation_observation"
TEMP = Path("<TMP_DIR>/prn09_comparable_v3_batches")
CAL_IQ = TEMP / "prn09_isolated_current_binary_145s_iq16.bin"
EXPECTED_IQ_BYTES = 14_499_900_000
CALIBRATION_SEED = 2026101100
CONFIRMATION_SEEDS = (2026101110, 2026101111, 2026101112,
                      2026101113, 2026101114)
SCENES = {
    "clean": (tuple(range(61, 77)), 68),
    "amp_0p05": (tuple(range(100, 116)), 108),
    "amp_0p10": (tuple(range(121, 137)), 129),
}
KEEP_EPOCH_RESPONSES = False


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run_command(command: list[str], log: Path, cwd: Path = ROOT) -> None:
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
        raise RuntimeError(f"refusing to overwrite IQ file: {output}")
    run_command([
        str(SIMULATOR), "--config", str(config), "--output", str(output),
        "--noise-seed", str(seed), "--multi-thread",
    ], log)
    if output.stat().st_size != EXPECTED_IQ_BYTES:
        raise RuntimeError(f"unexpected IQ size {output.stat().st_size}: {output}")
    if "Clipped samples: 0" not in log.read_text(errors="replace"):
        raise RuntimeError(f"simulator clipping check failed: {config}")


def prepare_state(run_dir: Path, iq: Path) -> Path:
    state_dir = run_dir / "state"
    state_dir.mkdir(parents=True, exist_ok=True)
    lines = STATE_TEMPLATE.read_text(encoding="utf-8").splitlines()
    replace_or_append(lines, "SignalSource.filename", str(iq))
    replace_or_append(lines, "Channel0.satellite", str(TARGET_PRN))
    replace_or_append(lines, "Tracking_1C.corr_scan_prn", str(TARGET_PRN))
    replace_or_append(
        lines, "Tracking_1C.state_dump_filename",
        str(state_dir / "tracking_state_prn09.csv"),
    )
    replace_or_append(lines, "Observables.dump", "false")
    replace_or_append(lines, "PVT.dump", "false")
    config = run_dir / "configs/state.conf"
    config.parent.mkdir(parents=True, exist_ok=True)
    config.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return config


def extract_scene(run_dir: Path, scene: str, windows: tuple[int, ...],
                  validation_iq: Path, state: Path) -> None:
    destination = run_dir / "open_loop" / scene
    destination.mkdir(parents=True, exist_ok=True)
    output = destination / "prn09_epoch_responses.csv"
    run_command([
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
    run_command([
        sys.executable, str(ROOT / "shared_ratio_existence_prn09_comparable.py"),
        "--run-dir", str(run_dir), "--banks-dir", str(BANKS),
        "--output", str(output),
    ], run_dir / "shared_ratio.log")
    return output


def apply_rule(records: list[dict], k: int, m: int) -> list[dict]:
    output = [dict(row) for row in records]
    for scene in SCENES:
        rows = sorted(
            (row for row in output if row["scene"] == scene),
            key=lambda row: int(row["window"]),
        )
        evidence = [bool(row["single_window_exists"]) for row in rows]
        for index, row in enumerate(rows):
            active = (
                index >= m - 1
                and sum(evidence[index - m + 1:index + 1]) >= k
            )
            row["gate_active"] = bool(active)
            row["rule_k"] = int(k)
            row["rule_m"] = int(m)
    return output


def summarize(records: list[dict], k: int, m: int) -> dict:
    gated = apply_rule(records, k, m)
    result = {"k": k, "m": m, "scenes": {}}
    for scene in SCENES:
        rows = [row for row in gated if row["scene"] == scene]
        active = sum(bool(row["gate_active"]) for row in rows)
        result["scenes"][scene] = {
            "windows": len(rows),
            "active": int(active),
            "active_rate": float(active / len(rows)),
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
        run_command([
            str(RECEIVER), f"--config_file={state_config}",
            "--logtostderr=1", "--colorlogtostderr=1",
        ], run_dir / "state_console.log", state_config.parent)
        state = run_dir / "state/tracking_state_prn09.csv"
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
            "seed": seed,
            "truth_accessed": False,
            "threshold": threshold,
            "rule": {"k": k, "m": m},
            "state_sha256": sha256(state),
            "scenes": summarize(records, k, m)["scenes"],
        }
        (run_dir / "persistence_summary.json").write_text(
            json.dumps(summary, indent=2) + "\n", encoding="utf-8"
        )
        if not KEEP_EPOCH_RESPONSES:
            for scene in SCENES:
                (run_dir / "open_loop" / scene / "prn09_epoch_responses.csv").unlink()
        return summary
    finally:
        if validation_iq.exists() and validation_iq.parent == TEMP:
            validation_iq.unlink()


def main(argv: list[str] | None = None) -> None:
    global KEEP_EPOCH_RESPONSES
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--keep-epoch-responses",
        action="store_true",
        help="retain 401-point mixed response CSVs for independent baselines",
    )
    args = parser.parse_args(argv)
    KEEP_EPOCH_RESPONSES = bool(args.keep_epoch_responses)
    for path in (
        SIMULATOR, RECEIVER, SCANNER, CAL_CONFIG, VALID_CONFIG, STATE_TEMPLATE,
        BANKS / "bank_68.npz", BANKS / "bank_108.npz", BANKS / "bank_129.npz",
        CAL_AUDIT, THRESHOLD_FILE, RULE_FILE, PROTOCOL,
    ):
        if not path.is_file():
            raise FileNotFoundError(path)
    audit = json.loads(CAL_AUDIT.read_text())
    if audit["target_prn"] != TARGET_PRN:
        raise RuntimeError("calibration target PRN mismatch")
    if audit["receiver_sha256"] != sha256(RECEIVER):
        raise RuntimeError("current receiver binary differs from PRN9 audit")
    threshold = float(json.loads(THRESHOLD_FILE.read_text())["T_null"])
    frozen = json.loads(RULE_FILE.read_text())
    selected = frozen["selection"]["selected"]
    if (int(selected["k"]), int(selected["m"])) != (3, 5):
        raise RuntimeError("frozen rule is not 3-of-5")
    k, m = 3, 5
    if CONFIRMATION_ROOT.exists() and any(CONFIRMATION_ROOT.iterdir()):
        raise RuntimeError(f"confirmation output is not empty: {CONFIRMATION_ROOT}")
    if TEMP.exists() and any(TEMP.iterdir()):
        raise RuntimeError(f"temporary directory is not empty: {TEMP}")

    TEMP.mkdir(parents=True, exist_ok=True)
    CONFIRMATION_ROOT.mkdir(parents=True, exist_ok=True)
    try:
        print("stage=generate_calibration_primary", flush=True)
        generate(CAL_CONFIG, CAL_IQ, CALIBRATION_SEED,
                 TEMP / "confirmation_calibration_primary_generation.log")
        if sha256(CAL_IQ) != audit["calibration_primary"]["sha256"]:
            raise RuntimeError("regenerated calibration IQ hash differs from audit")
        summaries = []
        for seed in CONFIRMATION_SEEDS:
            if shutil.disk_usage("/tmp").free < 18_000_000_000:
                raise RuntimeError("less than 18 GB available for next confirmation seed")
            summaries.append(process_seed(seed, threshold, k, m))
        aggregate = {}
        for scene in SCENES:
            active = sum(item["scenes"][scene]["active"] for item in summaries)
            windows = sum(item["scenes"][scene]["windows"] for item in summaries)
            aggregate[scene] = {
                "windows": windows,
                "active": active,
                "active_rate": float(active / windows),
            }
        payload = {
            "status": "prn09_comparable_shared_ratio_confirmation_complete",
            "truth_accessed": False,
            "target_prn": TARGET_PRN,
            "threshold": threshold,
            "threshold_file_sha256": sha256(THRESHOLD_FILE),
            "rule_file_sha256": sha256(RULE_FILE),
            "protocol_sha256": sha256(PROTOCOL),
            "calibration_audit_sha256": sha256(CAL_AUDIT),
            "rule": {"k": k, "m": m},
            "seeds": list(CONFIRMATION_SEEDS),
            "aggregate": aggregate,
            "per_seed": summaries,
        }
        (CONFIRMATION_ROOT / "CONFIRMATION_SUMMARY.json").write_text(
            json.dumps(payload, indent=2) + "\n", encoding="utf-8"
        )
        print(json.dumps({"rule": payload["rule"], "aggregate": aggregate},
                         indent=2), flush=True)
    finally:
        if CAL_IQ.exists() and CAL_IQ.parent == TEMP:
            CAL_IQ.unlink()
        if TEMP.exists() and not any(TEMP.iterdir()):
            TEMP.rmdir()


if __name__ == "__main__":
    main()
