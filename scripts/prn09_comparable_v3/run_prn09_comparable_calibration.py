#!/usr/bin/env python3
"""Calibrate PRN9 response banks with the current receiver binary."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

ROOT = Path(__file__).resolve().parent
PAPER_ROOT = ROOT.parent
sys.path.insert(0, str(PAPER_ROOT))

from receiver_chain_response_control import (  # noqa: E402
    ResponseSpec,
    build_bank,
    load_response_model,
)


TARGET_PRN = 9
WINDOWS = (68, 108, 129)
SPEC = ResponseSpec(prn=TARGET_PRN)
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
SIM_CONFIG = ROOT / "configs/signalsim_prn09_isolated_145s.json"
STATE_TEMPLATE = (
    PAPER_ROOT / "gps_l1_prn30_formal_pvt/final_sign_persistence/state/run.conf"
)
PLAN = ROOT / "CALIBRATION_PLAN.json"
CALIBRATION = ROOT / "calibration"
PRIMARY = Path(
    "<TMP_DIR>/prn09_comparable_v3_calibration/prn09_isolated_145s_iq16.bin"
)
STATE = CALIBRATION / "tracking_state_prn09.csv"
BANKS = ROOT / "response_banks"
AUDIT = ROOT / "CALIBRATION_AUDIT.json"
EXPECTED_IQ_BYTES = 14_499_900_000
CALIBRATION_SEED = 2026101100


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


def prepare_state_config() -> Path:
    CALIBRATION.mkdir(parents=True, exist_ok=True)
    lines = STATE_TEMPLATE.read_text(encoding="utf-8").splitlines()
    replace_or_append(lines, "SignalSource.filename", str(PRIMARY))
    replace_or_append(lines, "Channel0.satellite", str(TARGET_PRN))
    replace_or_append(lines, "Tracking_1C.corr_scan_prn", str(TARGET_PRN))
    replace_or_append(lines, "Tracking_1C.state_dump_filename", str(STATE))
    replace_or_append(lines, "Observables.dump", "false")
    replace_or_append(lines, "PVT.dump", "false")
    path = CALIBRATION / "state_prn09.conf"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def generate_primary() -> None:
    PRIMARY.parent.mkdir(parents=True, exist_ok=True)
    if PRIMARY.exists():
        raise RuntimeError(f"refusing to overwrite existing calibration IQ: {PRIMARY}")
    run_command([
        str(SIMULATOR), "--config", str(SIM_CONFIG), "--output", str(PRIMARY),
        "--noise-seed", str(CALIBRATION_SEED), "--multi-thread",
    ], CALIBRATION / "generation.log")
    if PRIMARY.stat().st_size != EXPECTED_IQ_BYTES:
        raise RuntimeError(f"unexpected calibration IQ size: {PRIMARY.stat().st_size}")
    if "Clipped samples: 0" not in (CALIBRATION / "generation.log").read_text(
        errors="replace"
    ):
        raise RuntimeError("calibration generation did not report zero clipping")


def run_state_receiver(config: Path) -> None:
    if STATE.exists():
        raise RuntimeError(f"refusing to overwrite tracking state: {STATE}")
    run_command([
        str(RECEIVER), f"--config_file={config}",
        "--logtostderr=1", "--colorlogtostderr=1",
    ], CALIBRATION / "state_console.log", CALIBRATION)
    if not STATE.is_file() or STATE.stat().st_size == 0:
        raise RuntimeError("PRN9 state receiver did not produce a state file")


def extract_gpu_reference(window: int) -> Path:
    destination = CALIBRATION / "gpu_reference" / f"window_{window}"
    destination.mkdir(parents=True, exist_ok=True)
    output = destination / "prn09_epoch_responses.csv"
    if output.exists():
        raise RuntimeError(f"refusing to overwrite GPU reference: {output}")
    run_command([
        str(SCANNER), "--prn", str(TARGET_PRN), "--state", str(STATE),
        "--primary", str(PRIMARY), "--secondary", str(PRIMARY),
        "--mixed", str(PRIMARY), "--output-dir", str(destination),
        "--epoch-output", str(output), "--epochs-per-window", "100",
        "--max-windows", "1", "--grid-ms", "1000",
        "--skip-ms", str(window * 1000), "--fs-hz", "25000000",
        "--if-hz", "-4580000",
    ], destination / "scan.log")
    rows = sum(1 for _ in output.open(encoding="utf-8"))
    if rows != 100 * 401 + 1:
        raise RuntimeError(f"window {window}: {rows} rows, expected {100 * 401 + 1}")
    return output


def gpu_shape_check(window: int, response_csv: Path) -> dict:
    by_epoch: dict[int, list[complex]] = {}
    with response_csv.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            epoch = int(row["epoch"])
            by_epoch.setdefault(epoch, []).append(
                complex(float(row["primary_i"]), float(row["primary_q"]))
            )
    if sorted(by_epoch) != list(range(100)):
        raise RuntimeError(f"window {window}: incomplete GPU epochs")
    reference = np.mean(
        np.asarray([by_epoch[epoch] for epoch in range(100)], complex), axis=0
    )
    if reference.shape != (401,):
        raise RuntimeError(f"window {window}: GPU reference shape {reference.shape}")
    model = load_response_model(BANKS / f"bank_{window}.npz", expected_spec=SPEC)
    predicted = model.exact([0.0])[0]
    gain = np.vdot(predicted, reference) / np.vdot(predicted, predicted)
    error = np.linalg.norm(reference - gain * predicted) / np.linalg.norm(reference)
    return {
        "window": window,
        "relative_shape_error": float(error),
        "gain_real": float(gain.real),
        "gain_imag": float(gain.imag),
        "pass_le_0p01": bool(error <= 0.01),
        "gpu_reference": str(response_csv),
        "gpu_reference_sha256": sha256(response_csv),
    }


def main() -> None:
    if AUDIT.exists():
        raise RuntimeError(f"calibration audit already exists: {AUDIT}")
    if BANKS.exists() and any(BANKS.iterdir()):
        raise RuntimeError(f"response bank directory is not empty: {BANKS}")
    for path in (SIMULATOR, RECEIVER, SCANNER, SIM_CONFIG, STATE_TEMPLATE, PLAN):
        if not path.is_file():
            raise FileNotFoundError(path)

    state_config = prepare_state_config()
    try:
        generate_primary()
        run_state_receiver(state_config)
        BANKS.mkdir(parents=True, exist_ok=True)
        checks = []
        bank_records = []
        for window in WINDOWS:
            bank_path = BANKS / f"bank_{window}.npz"
            build_bank(
                window, BANKS, manifest_path=PLAN, spec=SPEC,
                primary_path=PRIMARY, state_path=STATE,
                samples_per_window=25_000_000,
            )
            gpu_reference = extract_gpu_reference(window)
            check = gpu_shape_check(window, gpu_reference)
            checks.append(check)
            bank_records.append({
                "window": window,
                "path": str(bank_path),
                "bytes": bank_path.stat().st_size,
                "sha256": sha256(bank_path),
            })
            print(json.dumps(check), flush=True)
        if not all(check["pass_le_0p01"] for check in checks):
            raise RuntimeError("one or more PRN9 response banks failed GPU validation")

        payload = {
            "status": "prn09_comparable_calibration_complete_and_gpu_validated",
            "truth_accessed": False,
            "validation_iq_generated": False,
            "validation_outcomes_accessed": False,
            "target_prn": TARGET_PRN,
            "response_spec": SPEC.__dict__,
            "calibration_seed": CALIBRATION_SEED,
            "receiver": str(RECEIVER),
            "receiver_sha256": sha256(RECEIVER),
            "calibration_config": str(SIM_CONFIG),
            "calibration_config_sha256": sha256(SIM_CONFIG),
            "calibration_primary": {
                "path": str(PRIMARY),
                "bytes": PRIMARY.stat().st_size,
                "sha256": sha256(PRIMARY),
                "retained_after_audit": False,
            },
            "tracking_state": {
                "path": str(STATE),
                "bytes": STATE.stat().st_size,
                "sha256": sha256(STATE),
            },
            "state_config": {
                "path": str(state_config),
                "sha256": sha256(state_config),
            },
            "banks": bank_records,
            "gpu_forward_checks": checks,
            "inputs": {
                str(path): sha256(path)
                for path in (PLAN, SIM_CONFIG, STATE_TEMPLATE, RECEIVER, SCANNER)
            },
        }
        AUDIT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        print(json.dumps(payload, indent=2), flush=True)
    finally:
        if PRIMARY.exists() and PRIMARY.parent == Path("<TMP_DIR>/prn09_comparable_v3_calibration"):
            PRIMARY.unlink()


if __name__ == "__main__":
    main()
