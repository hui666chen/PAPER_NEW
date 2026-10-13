#!/usr/bin/env python3
"""Run the frozen PRN9 short-delay strong-reflection pseudorange pilot."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys
from zoneinfo import ZoneInfo

import numpy as np


WORKSPACE = Path("/home/zoxoai/Documents/ChatGPT/paper")
REPO = WORKSPACE / "PAPER_NEW"
RESULT_ROOT = REPO / "results/short_delay_pilot"
RAW_ROOT = WORKSPACE / "gps_l1_prn09_short_delay_pilot_20261013"
TEMP_ROOT = Path("/tmp/prn09_short_delay_pilot_20261013")

SIMULATOR = Path(
    "/home/zoxoai/chh/GPU/High/signalsim-new-multipath/build-if/IFdataGen"
)
RECEIVER = Path(
    "/home/zoxoai/chh/GPU/High/gnss-sdr-cuda-scan/"
    "build-cuda-scan/src/main/gnss-sdr"
)
SCANNER = Path(
    "/home/zoxoai/chh/GPU/High/gnss-sdr-cuda-scan/build-cuda-scan/"
    "src/algorithms/tracking/libs/frozen_reference_open_loop_scan"
)
STATE_TEMPLATE = (
    WORKSPACE / "gps_l1_prn30_formal_pvt/final_sign_persistence/state/run.conf"
)
RECEIVER_TEMPLATE = (
    WORKSPACE
    / "gps_l1_prn30_formal_pvt/final_sign_persistence/receiver/single/run.conf"
)
ASSISTANCE_DIR = (
    WORKSPACE / "gps_l1_prn30_formal_pvt/final_sign_persistence/receiver/single"
)

MIXED_CONFIG = RESULT_ROOT / "configs/signalsim_prn09_short_delay_mixed_96s_v2.json"
CLEAN_CONFIG = RESULT_ROOT / "configs/signalsim_prn09_short_delay_clean_96s_v2.json"
CAL_CONFIG = RESULT_ROOT / "configs/signalsim_prn09_short_delay_calibration_96s_v2.json"
PROTOCOL = RESULT_ROOT / "SHORT_DELAY_PILOT_PREREGISTRATION_20261013.md"
SEED_POLICY = RESULT_ROOT / "SEED_POLICY_20261013.md"
INFRASTRUCTURE_AMENDMENT = RESULT_ROOT / "INFRASTRUCTURE_AMENDMENT_20261013.md"
FREEZE_MANIFEST = RESULT_ROOT / "FROZEN_INPUTS_20261013.json"
COMPACT_RESULT = RESULT_ROOT / "SHORT_DELAY_PILOT_RESULTS_20261013.json"

sys.path.insert(0, str(WORKSPACE))
sys.path.insert(0, str(WORKSPACE / "gps_l1_prn09_shared_ratio_v2"))
sys.path.insert(0, str(REPO / "scripts/baselines"))

import shared_ratio_existence as dual_base  # noqa: E402
from gps_l1_prn30_formal_pvt.build_common_observable_rinex import (  # noqa: E402
    parse_epoch,
    relative_second,
)
from gps_l1_prn30_formal_pvt.estimate_formal_windows import (  # noqa: E402
    load_windows,
)
from monte_carlo_resolution_baseline import build_dictionary  # noqa: E402
from receiver_chain_response_control import (  # noqa: E402
    ResponseSpec,
    build_bank,
    load_response_model,
)
from correlator_baselines import (  # noqa: E402
    BaselineConfig,
    double_delta_estimate,
    eml_estimate,
    read_epoch_scans,
)


TARGET_PRN = 9
TARGET_SATELLITE = "G09"
SPEC = ResponseSpec(prn=TARGET_PRN)
SAMPLE_RATE_HZ = 25_000_000
EXPECTED_IQ_BYTES = 9_599_900_000
CHIP_METRES = 299_792_458.0 / 1.023e6
CALIBRATION_SEED = 2026101300
VALIDATION_SEEDS = (2026101311, 2026101312, 2026101313, 2026101314, 2026101320)
BOOTSTRAP_SEED = 2026101399
BOOTSTRAP_REPLICATES = 200_000
METHODS = ("raw", "eml", "double_delta", "dual")
CELLS = {
    "d005_a025": {
        "windows": tuple(range(36, 44)),
        "bank_window": 39,
        "delay_chip": 0.05,
        "amplitude_ratio": 0.25,
        "segment": (34.0, 46.0),
    },
    "d005_a050": {
        "windows": tuple(range(50, 58)),
        "bank_window": 53,
        "delay_chip": 0.05,
        "amplitude_ratio": 0.50,
        "segment": (48.0, 60.0),
    },
    "d010_a025": {
        "windows": tuple(range(64, 72)),
        "bank_window": 67,
        "delay_chip": 0.10,
        "amplitude_ratio": 0.25,
        "segment": (62.0, 74.0),
    },
    "d010_a050": {
        "windows": tuple(range(78, 86)),
        "bank_window": 81,
        "delay_chip": 0.10,
        "amplitude_ratio": 0.50,
        "segment": (76.0, 88.0),
    },
}
BASELINE_CONFIG = BaselineConfig()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run_command(command: list[str], log: Path, cwd: Path) -> None:
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("w", encoding="utf-8", errors="replace") as stream:
        subprocess.run(
            command,
            cwd=cwd,
            stdout=stream,
            stderr=subprocess.STDOUT,
            check=True,
        )


def replace_or_append(lines: list[str], key: str, value: str) -> None:
    prefix = key + "="
    for index, line in enumerate(lines):
        if line.startswith(prefix):
            lines[index] = prefix + value
            return
    lines.append(prefix + value)


def verify_frozen_inputs() -> dict:
    frozen = json.loads(FREEZE_MANIFEST.read_text(encoding="utf-8"))
    if frozen["status"] != "frozen_before_validation_generation":
        raise RuntimeError("invalid freeze status")
    for path_text, expected in frozen["sha256"].items():
        path = Path(path_text)
        actual = sha256(path)
        if actual != expected:
            raise RuntimeError(f"frozen input changed: {path}\n{actual}\n{expected}")
    return frozen


def generate(config: Path, output: Path, seed: int, log: Path) -> None:
    if output.exists():
        raise RuntimeError(f"refusing to overwrite IQ: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    run_command(
        [
            str(SIMULATOR),
            "--config",
            str(config),
            "--output",
            str(output),
            "--noise-seed",
            str(seed),
            "--multi-thread",
        ],
        log,
        RAW_ROOT,
    )
    if output.stat().st_size != EXPECTED_IQ_BYTES:
        raise RuntimeError(
            f"unexpected IQ bytes {output.stat().st_size}, expected {EXPECTED_IQ_BYTES}"
        )
    text = log.read_text(encoding="utf-8", errors="replace")
    if "Clipped samples: 0" not in text:
        raise RuntimeError(f"clipping or missing clipping record: {log}")


def prepare_state_config(directory: Path, iq: Path, label: str) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    state = directory / f"tracking_state_prn09_{label}.csv"
    lines = STATE_TEMPLATE.read_text(encoding="utf-8").splitlines()
    replace_or_append(lines, "SignalSource.filename", str(iq))
    replace_or_append(lines, "SignalSource.samples", "0")
    replace_or_append(lines, "Channel0.satellite", str(TARGET_PRN))
    replace_or_append(lines, "Tracking_1C.corr_scan_prn", str(TARGET_PRN))
    replace_or_append(lines, "Tracking_1C.state_dump_filename", str(state))
    replace_or_append(lines, "Observables.dump", "false")
    replace_or_append(lines, "PVT.dump", "false")
    config = directory / f"state_{label}.conf"
    config.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return config


def run_state_receiver(directory: Path, iq: Path, label: str) -> Path:
    config = prepare_state_config(directory, iq, label)
    state = directory / f"tracking_state_prn09_{label}.csv"
    run_command(
        [
            str(RECEIVER),
            f"--config_file={config}",
            "--logtostderr=1",
            "--colorlogtostderr=1",
        ],
        directory / f"state_{label}.log",
        directory,
    )
    if not state.is_file() or state.stat().st_size == 0:
        raise RuntimeError(f"missing tracking state: {state}")
    return state


def write_zero_schedule(path: Path) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=("start_sample", "end_sample", "prn", "correction_m"),
        )
        writer.writeheader()
        for second in range(96):
            writer.writerow(
                {
                    "start_sample": second * SAMPLE_RATE_HZ,
                    "end_sample": (second + 1) * SAMPLE_RATE_HZ - 1,
                    "prn": TARGET_PRN,
                    "correction_m": "0.000000000000",
                }
            )


def prepare_observable_receiver(directory: Path, iq: Path) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    zero_schedule = directory / "zero_corrections.csv"
    write_zero_schedule(zero_schedule)
    lines = RECEIVER_TEMPLATE.read_text(encoding="utf-8").splitlines()
    replace_or_append(lines, "SignalSource.filename", str(iq))
    replace_or_append(lines, "SignalSource.samples", "0")
    replace_or_append(lines, "GNSS-SDR.SUPL_gps_enabled", "true")
    replace_or_append(lines, "GNSS-SDR.SUPL_read_gps_assistance_xml", "true")
    replace_or_append(
        lines,
        "GNSS-SDR.SUPL_gps_ephemeris_xml",
        str(ASSISTANCE_DIR / "gps_ephemeris.xml"),
    )
    replace_or_append(
        lines, "GNSS-SDR.SUPL_gps_iono_xml", str(ASSISTANCE_DIR / "gps_iono.xml")
    )
    replace_or_append(
        lines,
        "GNSS-SDR.SUPL_gps_utc_model_xml",
        str(ASSISTANCE_DIR / "gps_utc_model.xml"),
    )
    replace_or_append(lines, "Observables.dump_filename", str(directory / "observables"))
    replace_or_append(lines, "PVT.dump_filename", str(directory / "pvt"))
    replace_or_append(lines, "PVT.output_path", str(directory))
    replace_or_append(lines, "Observables.pseudorange_correction_file", str(zero_schedule))
    replace_or_append(lines, "PVT.output_rate_ms", "100")
    replace_or_append(lines, "PVT.display_rate_ms", "500")
    config = directory / "receiver.conf"
    config.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return config


def run_observable_receiver(directory: Path, iq: Path) -> tuple[Path, Path]:
    config = prepare_observable_receiver(directory, iq)
    run_command(
        [
            str(RECEIVER),
            f"--config_file={config}",
            "--logtostderr=1",
            "--colorlogtostderr=1",
        ],
        directory / "console.log",
        directory,
    )
    observations = sorted(directory.glob("*.??O"))
    navigation = sorted(directory.glob("*.??N"))
    if len(observations) != 1 or len(navigation) != 1:
        raise RuntimeError(
            f"expected one observation/navigation RINEX in {directory}: "
            f"{observations}, {navigation}"
        )
    return observations[0], navigation[0]


def extract_scene(
    root: Path,
    scene: str,
    windows: tuple[int, ...],
    primary_iq: Path,
    mixed_iq: Path,
    state: Path,
) -> tuple[Path, Path]:
    destination = root / "open_loop" / scene
    destination.mkdir(parents=True, exist_ok=True)
    epoch_output = destination / "prn09_epoch_responses.csv"
    run_command(
        [
            str(SCANNER),
            "--prn",
            str(TARGET_PRN),
            "--state",
            str(state),
            "--primary",
            str(primary_iq),
            "--secondary",
            str(primary_iq),
            "--mixed",
            str(mixed_iq),
            "--output-dir",
            str(destination),
            "--epoch-output",
            str(epoch_output),
            "--epochs-per-window",
            "100",
            "--max-windows",
            str(len(windows)),
            "--grid-ms",
            "1000",
            "--skip-ms",
            str(windows[0] * 1000),
            "--fs-hz",
            "25000000",
            "--if-hz",
            "-4580000",
        ],
        destination / "scan.log",
        destination,
    )
    expected_epoch_rows = len(windows) * 100 * 401 + 1
    if sum(1 for _ in epoch_output.open(encoding="utf-8")) != expected_epoch_rows:
        raise RuntimeError(f"incomplete epoch response: {epoch_output}")
    average_output = destination / "frozen_open_loop_prn9_responses.csv"
    expected_average_rows = len(windows) * 401 + 1
    if not average_output.is_file() or sum(
        1 for _ in average_output.open(encoding="utf-8")
    ) != expected_average_rows:
        raise RuntimeError(f"incomplete average response: {average_output}")
    return epoch_output, average_output


def baseline_estimates(path: Path) -> dict[str, dict[int, dict]]:
    scans = read_epoch_scans(path)
    output: dict[str, dict[int, dict]] = {"eml": {}, "double_delta": {}}
    for scan in scans:
        output["eml"][scan.window] = eml_estimate(scan, BASELINE_CONFIG)
        output["double_delta"][scan.window] = double_delta_estimate(
            scan, BASELINE_CONFIG
        )
    return output


def dual_estimates(
    epoch_path: Path, windows: tuple[int, ...], bank_path: Path
) -> dict[int, dict]:
    observations = load_windows(epoch_path, windows)
    model = load_response_model(bank_path, expected_spec=SPEC)
    delay_grid, dictionary, _, _, _ = build_dictionary(model, SPEC)
    rho_values = dual_base.make_rho_grid()
    output = {}
    for window in windows:
        values, weights, magnitude = observations[window]
        record = dual_base.analyze_window(
            values, weights, dictionary, delay_grid, rho_values
        )
        fits = (record["A_to_B"]["h1"], record["B_to_A"]["h1"])
        record["tau_chip"] = float(np.mean([row["tau_chip"] for row in fits]))
        record["tau_boundary"] = bool(
            any(
                np.isclose(row["tau_chip"], dual_base.TAU[0])
                or np.isclose(row["tau_chip"], dual_base.TAU[-1])
                for row in fits
            )
        )
        record["separation_boundary"] = bool(
            any(
                np.isclose(row["separation_chip"], dual_base.SEPARATION[0])
                or np.isclose(row["separation_chip"], dual_base.SEPARATION[-1])
                for row in fits
            )
        )
        record["prompt_magnitude_median"] = float(np.median(magnitude))
        output[window] = record
    return output


def calibration_stage() -> dict:
    directory = RAW_ROOT / "calibration_v2_after_acquisition_amendment"
    manifest_path = directory / "CALIBRATION_RESULT.json"
    cal_iq = TEMP_ROOT / "calibration_seed_2026101300_96s_v2_iq16.bin"
    if manifest_path.is_file():
        if not cal_iq.is_file() or cal_iq.stat().st_size != EXPECTED_IQ_BYTES:
            raise RuntimeError(
                "calibration result exists but its temporary IQ is missing; "
                "remove the incomplete raw calibration directory before restarting"
            )
        return json.loads(manifest_path.read_text(encoding="utf-8"))
    generation_log = directory / "generation.log"
    if directory.exists() and any(directory.iterdir()):
        resumable = (
            cal_iq.is_file()
            and cal_iq.stat().st_size == EXPECTED_IQ_BYTES
            and generation_log.is_file()
            and "Clipped samples: 0"
            in generation_log.read_text(encoding="utf-8", errors="replace")
        )
        if not resumable:
            raise RuntimeError(f"incomplete non-resumable calibration directory: {directory}")
        print("stage=calibration_resume_after_size_check_erratum", flush=True)
    else:
        directory.mkdir(parents=True, exist_ok=True)
        print("stage=calibration_generate", flush=True)
        generate(CAL_CONFIG, cal_iq, CALIBRATION_SEED, generation_log)
    state = run_state_receiver(directory / "state", cal_iq, "calibration")
    banks = directory / "response_banks"
    responses = {}
    for cell, design in CELLS.items():
        print(f"stage=calibration_scan cell={cell}", flush=True)
        epoch, average = extract_scene(
            directory, cell, design["windows"], cal_iq, cal_iq, state
        )
        build_bank(
            design["bank_window"],
            banks,
            manifest_path=CAL_CONFIG,
            spec=SPEC,
            primary_path=cal_iq,
            state_path=state,
            samples_per_window=SAMPLE_RATE_HZ,
        )
        responses[cell] = {"epoch": epoch, "average": average}

    offsets = {}
    calibration_records = {}
    for cell, design in CELLS.items():
        base = baseline_estimates(responses[cell]["average"])
        dual = dual_estimates(
            responses[cell]["epoch"],
            design["windows"],
            banks / f"bank_{design['bank_window']}.npz",
        )
        cell_offsets = {}
        cell_records = {"eml": [], "double_delta": [], "dual": []}
        for method in ("eml", "double_delta"):
            valid = []
            for window in design["windows"]:
                row = base[method][window]
                item = {
                    "window": window,
                    "tau_chip": float(row["tau_chip"]),
                    "boundary": bool(row["boundary"]),
                }
                cell_records[method].append(item)
                if not item["boundary"]:
                    valid.append(item["tau_chip"])
            if len(valid) < 4:
                raise RuntimeError(f"{cell}/{method}: fewer than 4 valid calibration windows")
            cell_offsets[method] = float(np.median(valid))
        valid_dual = []
        for window in design["windows"]:
            row = dual[window]
            item = {
                "window": window,
                "tau_chip": float(row["tau_chip"]),
                "tau_boundary": bool(row["tau_boundary"]),
                "separation_boundary": bool(row["separation_boundary"]),
            }
            cell_records["dual"].append(item)
            if not item["tau_boundary"]:
                valid_dual.append(item["tau_chip"])
        if len(valid_dual) < 4:
            raise RuntimeError(f"{cell}/dual: fewer than 4 valid calibration windows")
        cell_offsets["dual"] = float(np.median(valid_dual))
        offsets[cell] = cell_offsets
        calibration_records[cell] = cell_records

    result = {
        "status": "independent_calibration_complete",
        "seed": CALIBRATION_SEED,
        "truth_accessed": False,
        "iq_sha256": sha256(cal_iq),
        "state_sha256": sha256(state),
        "offsets_chip": offsets,
        "records": calibration_records,
        "banks": {
            cell: {
                "window": design["bank_window"],
                "sha256": sha256(banks / f"bank_{design['bank_window']}.npz"),
            }
            for cell, design in CELLS.items()
        },
    }
    manifest_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result


def iq_block_rms(mixed: np.memmap, clean: np.memmap, second: float) -> float:
    count = 200_000
    center = int(round(second * SAMPLE_RATE_HZ))
    first = center - count // 2
    delta = mixed[first : first + count].astype(np.float64) - clean[
        first : first + count
    ].astype(np.float64)
    return float(np.sqrt(np.mean(delta * delta)))


def injection_qc(mixed_iq: Path, clean_iq: Path, mixed_log: Path, clean_log: Path) -> dict:
    mixed_text = mixed_log.read_text(encoding="utf-8", errors="replace")
    clean_text = clean_log.read_text(encoding="utf-8", errors="replace")
    if "Loaded 4 time-varying multipath segment(s)" not in mixed_text:
        raise RuntimeError("mixed simulator log did not load all four multipath segments")
    if "multipath segment(s)" in clean_text:
        raise RuntimeError("counterfactual clean simulator log loaded multipath")
    mixed = np.memmap(mixed_iq, dtype="<i2", mode="r").reshape(-1, 2)
    clean = np.memmap(clean_iq, dtype="<i2", mode="r").reshape(-1, 2)
    guard = [iq_block_rms(mixed, clean, value) for value in (5.0, 93.0)]
    guard_max = max(guard)
    cells = {}
    for cell, design in CELLS.items():
        start, end = design["segment"]
        interior = iq_block_rms(mixed, clean, 0.5 * (start + end))
        ratio = interior / max(guard_max, 1.0e-30)
        passed = interior > 0.0 and ratio >= 1000.0
        cells[cell] = {
            "interior_difference_rms_iq16": interior,
            "guard_max_difference_rms_iq16": guard_max,
            "interior_to_guard_ratio": ratio,
            "pass": passed,
        }
        if not passed:
            raise RuntimeError(f"{cell}: paired-IQ injection hard gate failed")
    del mixed, clean
    return {
        "status": "PASS",
        "mixed_loaded_four_segments": True,
        "clean_loaded_no_segments": True,
        "guard_rms_iq16": guard,
        "cells": cells,
    }


def extract_rinex_c1c(path: Path) -> dict[int, dict[str, float]]:
    values: dict[int, dict[str, float]] = {}
    second = None
    in_header = True
    with path.open(encoding="ascii") as stream:
        for line in stream:
            if in_header:
                if "END OF HEADER" in line:
                    in_header = False
                continue
            if line.startswith(">"):
                second = relative_second(parse_epoch(line))
                values.setdefault(second, {})
                continue
            if second is None or not line.startswith("G"):
                continue
            text = line[3:17]
            if text.strip():
                values[second][line[:3]] = float(text)
    return values


def counterfactual_errors(
    mixed_path: Path, clean_path: Path, windows: tuple[int, ...]
) -> dict[int, dict]:
    mixed = extract_rinex_c1c(mixed_path)
    clean = extract_rinex_c1c(clean_path)
    output = {}
    for second in windows:
        if TARGET_SATELLITE not in mixed.get(second, {}):
            raise RuntimeError(f"mixed RINEX missing {TARGET_SATELLITE} at {second}")
        if TARGET_SATELLITE not in clean.get(second, {}):
            raise RuntimeError(f"clean RINEX missing {TARGET_SATELLITE} at {second}")
        common = sorted((set(mixed[second]) & set(clean[second])) - {TARGET_SATELLITE})
        if len(common) < 4:
            raise RuntimeError(f"fewer than four common non-target satellites at {second}")
        common_differences = np.asarray(
            [mixed[second][sat] - clean[second][sat] for sat in common], dtype=float
        )
        common_clock = float(np.median(common_differences))
        target_difference = float(
            mixed[second][TARGET_SATELLITE] - clean[second][TARGET_SATELLITE]
        )
        output[second] = {
            "raw_error_m": target_difference - common_clock,
            "target_difference_m": target_difference,
            "common_clock_difference_m": common_clock,
            "common_non_target_satellites": common,
        }
    return output


def method_metric(errors: list[float]) -> dict:
    values = np.asarray(errors, dtype=float)
    return {
        "seconds": int(len(values)),
        "rmse_m": float(np.sqrt(np.mean(values * values))),
        "mae_m": float(np.mean(np.abs(values))),
        "median_error_m": float(np.median(values)),
        "p95_absolute_error_m": float(np.quantile(np.abs(values), 0.95)),
        "gross_error_rate_gt_10m": float(np.mean(np.abs(values) > 10.0)),
    }


def process_seed(seed: int, calibration: dict) -> dict:
    run_dir = RAW_ROOT / f"seed_{seed}"
    result_path = run_dir / "SEED_RESULT.json"
    if result_path.is_file():
        return json.loads(result_path.read_text(encoding="utf-8"))
    if run_dir.exists() and any(run_dir.iterdir()):
        raise RuntimeError(f"incomplete non-empty run directory: {run_dir}")
    run_dir.mkdir(parents=True, exist_ok=True)
    mixed_iq = TEMP_ROOT / f"mixed_seed_{seed}_iq16.bin"
    clean_iq = TEMP_ROOT / f"clean_seed_{seed}_iq16.bin"
    cal_iq = TEMP_ROOT / "calibration_seed_2026101300_96s_v2_iq16.bin"
    if not cal_iq.is_file():
        raise RuntimeError("independent calibration IQ is missing")
    try:
        print(f"seed={seed} stage=generate_mixed", flush=True)
        generate(MIXED_CONFIG, mixed_iq, seed, run_dir / "mixed_generation.log")
        print(f"seed={seed} stage=generate_clean", flush=True)
        generate(CLEAN_CONFIG, clean_iq, seed, run_dir / "clean_generation.log")
        qc = injection_qc(
            mixed_iq,
            clean_iq,
            run_dir / "mixed_generation.log",
            run_dir / "clean_generation.log",
        )
        print(f"seed={seed} stage=mixed_tracking_state", flush=True)
        mixed_state = run_state_receiver(run_dir / "mixed_state", mixed_iq, "mixed")

        estimated = {}
        for cell, design in CELLS.items():
            print(f"seed={seed} stage=response cell={cell}", flush=True)
            epoch, average = extract_scene(
                run_dir,
                cell,
                design["windows"],
                cal_iq,
                mixed_iq,
                mixed_state,
            )
            base = baseline_estimates(average)
            dual = dual_estimates(
                epoch,
                design["windows"],
                RAW_ROOT
                / "calibration_v2_after_acquisition_amendment/response_banks"
                / f"bank_{design['bank_window']}.npz",
            )
            estimated[cell] = {
                "eml": base["eml"],
                "double_delta": base["double_delta"],
                "dual": dual,
            }

        print(f"seed={seed} stage=receiver_mixed", flush=True)
        mixed_obs, mixed_nav = run_observable_receiver(
            run_dir / "receiver_mixed", mixed_iq
        )
        print(f"seed={seed} stage=receiver_clean", flush=True)
        clean_obs, clean_nav = run_observable_receiver(
            run_dir / "receiver_clean", clean_iq
        )

        rows = []
        metrics = {}
        for cell, design in CELLS.items():
            truth = counterfactual_errors(mixed_obs, clean_obs, design["windows"])
            by_method = {method: [] for method in METHODS}
            boundary_counts = {method: 0 for method in METHODS}
            for second in design["windows"]:
                raw_error = float(truth[second]["raw_error_m"])
                by_method["raw"].append(raw_error)
                rows.append(
                    {
                        "seed": seed,
                        "cell": cell,
                        "second": second,
                        "method": "raw",
                        "tau_hat_chip": None,
                        "zero_chip": 0.0,
                        "correction_m": 0.0,
                        "raw_error_m": raw_error,
                        "residual_error_m": raw_error,
                        "boundary": False,
                        **truth[second],
                    }
                )
                for method in ("eml", "double_delta", "dual"):
                    estimate = estimated[cell][method][second]
                    tau = float(estimate["tau_chip"])
                    zero = float(calibration["offsets_chip"][cell][method])
                    correction = -(tau - zero) * CHIP_METRES
                    residual = raw_error + correction
                    boundary = bool(
                        estimate["boundary"]
                        if method in ("eml", "double_delta")
                        else estimate["tau_boundary"] or estimate["separation_boundary"]
                    )
                    boundary_counts[method] += int(boundary)
                    by_method[method].append(residual)
                    row = {
                        "seed": seed,
                        "cell": cell,
                        "second": second,
                        "method": method,
                        "tau_hat_chip": tau,
                        "zero_chip": zero,
                        "correction_m": correction,
                        "raw_error_m": raw_error,
                        "residual_error_m": residual,
                        "boundary": boundary,
                        **truth[second],
                    }
                    if method == "dual":
                        row.update(
                            {
                                "separation_a_chip": estimate["A_to_B"]["h1"][
                                    "separation_chip"
                                ],
                                "separation_b_chip": estimate["B_to_A"]["h1"][
                                    "separation_chip"
                                ],
                                "rho_a": estimate["A_to_B"]["h1"]["rho_magnitude"],
                                "rho_b": estimate["B_to_A"]["h1"]["rho_magnitude"],
                            }
                        )
                    rows.append(row)
            metrics[cell] = {}
            for method in METHODS:
                item = method_metric(by_method[method])
                item["boundary_rate"] = boundary_counts[method] / len(design["windows"])
                metrics[cell][method] = item

        with (run_dir / "epoch_scores.csv").open(
            "w", newline="", encoding="utf-8"
        ) as stream:
            fields = sorted({key for row in rows for key in row})
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)
        result = {
            "status": "seed_complete",
            "seed": seed,
            "truth_used_for_estimation": False,
            "counterfactual_clean_used_for_scoring_only": True,
            "injection_qc": qc,
            "mixed_iq_sha256": sha256(mixed_iq),
            "clean_iq_sha256": sha256(clean_iq),
            "mixed_state_sha256": sha256(mixed_state),
            "mixed_observation_sha256": sha256(mixed_obs),
            "clean_observation_sha256": sha256(clean_obs),
            "mixed_navigation_sha256": sha256(mixed_nav),
            "clean_navigation_sha256": sha256(clean_nav),
            "metrics": metrics,
            "epoch_score_sha256": sha256(run_dir / "epoch_scores.csv"),
        }
        result_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
        return result
    finally:
        for path in (mixed_iq, clean_iq):
            if path.exists() and path.parent == TEMP_ROOT:
                path.unlink()


def paired_interval(values: np.ndarray, indices: np.ndarray) -> list[float]:
    means = np.mean(values[indices], axis=1)
    return [float(value) for value in np.quantile(means, [0.025, 0.975])]


def aggregate(results: list[dict], calibration: dict, frozen: dict) -> dict:
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    indices = rng.integers(
        0, len(results), size=(BOOTSTRAP_REPLICATES, len(results))
    )
    cells = {}
    any_pass = False
    for cell in CELLS:
        method_summary = {}
        for method in METHODS:
            values = np.asarray(
                [result["metrics"][cell][method]["rmse_m"] for result in results],
                dtype=float,
            )
            method_summary[method] = {
                "mean_seed_rmse_m": float(np.mean(values)),
                "per_seed_rmse_m": values.tolist(),
            }
        dd = np.asarray(method_summary["double_delta"]["per_seed_rmse_m"])
        dual = np.asarray(method_summary["dual"]["per_seed_rmse_m"])
        absolute = dd - dual
        relative = 1.0 - dual / dd
        mean_reduction = 1.0 - float(np.mean(dual)) / float(np.mean(dd))
        absolute_ci = paired_interval(absolute, indices)
        relative_ci = paired_interval(relative, indices)
        passed = bool(mean_reduction >= 0.30 and absolute_ci[0] > 0.0)
        any_pass = any_pass or passed
        cells[cell] = {
            "methods": method_summary,
            "dual_vs_dd": {
                "mean_reduction_fraction_from_equal_seed_rmse_means": mean_reduction,
                "mean_paired_absolute_benefit_m": float(np.mean(absolute)),
                "absolute_benefit_ci95_m": absolute_ci,
                "mean_paired_relative_benefit": float(np.mean(relative)),
                "relative_benefit_ci95": relative_ci,
                "pass": passed,
            },
        }
    return {
        "status": "short_delay_pilot_complete",
        "completed_at_local": datetime.now(ZoneInfo("Asia/Shanghai")).isoformat(),
        "decision": "PASS_EXPAND" if any_pass else "FAIL_STOP_SUPERIORITY_CLAIM",
        "pass_rule_satisfied_by_at_least_one_cell": any_pass,
        "seeds": list(VALIDATION_SEEDS),
        "bootstrap_seed": BOOTSTRAP_SEED,
        "bootstrap_replicates": BOOTSTRAP_REPLICATES,
        "truth_used_for_estimation": False,
        "counterfactual_clean_used_for_scoring_only": True,
        "frozen_input_manifest_sha256": sha256(FREEZE_MANIFEST),
        "calibration": calibration,
        "cells": cells,
        "per_seed": results,
    }


def parse_seeds(text: str) -> tuple[int, ...]:
    values = tuple(int(value) for value in text.split(",") if value.strip())
    if not values or any(value not in VALIDATION_SEEDS for value in values):
        raise ValueError(f"seeds must be a subset of {VALIDATION_SEEDS}")
    return values


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--seeds",
        default=",".join(str(seed) for seed in VALIDATION_SEEDS),
        help="comma-separated frozen validation seeds; completed seeds are reused",
    )
    parser.add_argument(
        "--keep-calibration-iq",
        action="store_true",
        help="retain the reproducible 9.6 GB calibration IQ after all seeds complete",
    )
    args = parser.parse_args(argv)
    selected_seeds = parse_seeds(args.seeds)
    required = (
        SIMULATOR,
        RECEIVER,
        SCANNER,
        STATE_TEMPLATE,
        RECEIVER_TEMPLATE,
        MIXED_CONFIG,
        CLEAN_CONFIG,
        CAL_CONFIG,
        PROTOCOL,
        SEED_POLICY,
        INFRASTRUCTURE_AMENDMENT,
        FREEZE_MANIFEST,
    )
    for path in required:
        if not path.is_file():
            raise FileNotFoundError(path)
    frozen = verify_frozen_inputs()
    RAW_ROOT.mkdir(parents=True, exist_ok=True)
    TEMP_ROOT.mkdir(parents=True, exist_ok=True)
    if shutil.disk_usage(TEMP_ROOT).free < 30_000_000_000:
        raise RuntimeError("fewer than 30 GB free before short-delay pilot")
    calibration = calibration_stage()
    for seed in selected_seeds:
        if shutil.disk_usage(TEMP_ROOT).free < 18_000_000_000:
            raise RuntimeError("fewer than 18 GB free before next paired validation seed")
        process_seed(seed, calibration)

    completed = []
    for seed in VALIDATION_SEEDS:
        path = RAW_ROOT / f"seed_{seed}/SEED_RESULT.json"
        if path.is_file():
            completed.append(json.loads(path.read_text(encoding="utf-8")))
    if len(completed) == len(VALIDATION_SEEDS):
        result = aggregate(completed, calibration, frozen)
        COMPACT_RESULT.write_text(
            json.dumps(result, indent=2) + "\n", encoding="utf-8"
        )
        if not args.keep_calibration_iq:
            cal_iq = TEMP_ROOT / "calibration_seed_2026101300_96s_v2_iq16.bin"
            if cal_iq.exists() and cal_iq.parent == TEMP_ROOT:
                cal_iq.unlink()
        print(
            json.dumps(
                {
                    "decision": result["decision"],
                    "cells": {
                        cell: row["dual_vs_dd"] for cell, row in result["cells"].items()
                    },
                },
                indent=2,
            ),
            flush=True,
        )
    else:
        print(
            json.dumps(
                {
                    "status": "partial_complete",
                    "completed_seeds": [row["seed"] for row in completed],
                    "remaining_seeds": [
                        seed
                        for seed in VALIDATION_SEEDS
                        if seed not in {row["seed"] for row in completed}
                    ],
                },
                indent=2,
            ),
            flush=True,
        )


if __name__ == "__main__":
    main()
