#!/usr/bin/env python3
"""Run the frozen short PRN28 tracking precheck."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parent
PAPER_ROOT = ROOT.parent
PROTOCOL = ROOT / "PRN28_MID_ELEVATION_PROTOCOL.md"
CONFIG = ROOT / "configs/signalsim_prn28_precheck_20s.json"
SIGNALSIM = Path("<LOCAL_USER>/chh/GNSS/SignalSim-main")
RINEX = SIGNALSIM / "EphData/BRDC00IGS_R_20211700000_01D_MN.rnx"
SIMULATOR = Path("<LOCAL_USER>/chh/GPU/High/signalsim-new-multipath/build-if/IFdataGen")
RECEIVER = Path("<LOCAL_USER>/chh/GPU/High/gnss-sdr-cuda-scan/build-cuda-scan/src/main/gnss-sdr")
STATE_TEMPLATE = PAPER_ROOT / "gps_l1_prn30_formal_pvt/final_sign_persistence/state/run.conf"
SOURCE_GEOMETRY = PAPER_ROOT / "extract_signalsim_gps_geometry_at_time.cpp"
OUT = ROOT / "precheck"
TMP = Path("<TMP_DIR>/prn28_mid_elevation_precheck")
TARGET_PRN = 28
SEED = 2026101117
EXPECTED_IQ_BYTES = 1_999_900_000
FROZEN_EPOCH = (2021, 6, 19, 12, 25, 40)
EXPECTED_VISIBLE = (5, 7, 8, 9, 13, 14, 17, 20, 21, 28, 30)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run(command: list[str], log: Path, cwd: Path = ROOT) -> None:
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("w", encoding="utf-8", errors="replace") as stream:
        subprocess.run(command, cwd=cwd, stdout=stream,
                       stderr=subprocess.STDOUT, check=True)


def extract_geometry() -> list[dict]:
    binary = Path(tempfile.gettempdir()) / "extract_signalsim_gps_geometry_at_time_prn28"
    sources = ("NavData.cpp", "Rinex.cpp", "GnssTime.cpp", "Coordinate.cpp",
               "Almanac.cpp", "FastMath.cpp", "MessageOutput.cpp")
    command = ["g++", "-std=c++11", "-O2", f"-I{SIGNALSIM / 'inc'}",
               str(SOURCE_GEOMETRY)]
    command.extend(str(SIGNALSIM / "src" / name) for name in sources)
    command.extend(["-o", str(binary)])
    subprocess.run(command, check=True, capture_output=True, text=True)
    result = subprocess.run(
        [str(binary), str(RINEX), *map(str, FROZEN_EPOCH)],
        check=True, capture_output=True, text=True,
    )
    rows = [{
        "prn": int(row["prn"]),
        "elevation_deg": float(row["elevation_deg"]),
        "azimuth_deg": float(row["azimuth_deg"]),
        "los_ecef": [float(row[key]) for key in ("los_x", "los_y", "los_z")],
    } for row in csv.DictReader(result.stdout.splitlines())]
    if tuple(row["prn"] for row in rows) != EXPECTED_VISIBLE:
        raise RuntimeError("frozen visible-satellite set changed")
    target = next(row for row in rows if row["prn"] == TARGET_PRN)
    if not 40.0 <= target["elevation_deg"] <= 50.0:
        raise RuntimeError(f"target outside 40-50 degrees: {target}")
    return rows


def prepare_state(iq: Path, state: Path) -> Path:
    lines = STATE_TEMPLATE.read_text(encoding="utf-8").splitlines()
    values = {
        "SignalSource.filename": str(iq),
        "Channel0.satellite": str(TARGET_PRN),
        "Tracking_1C.corr_scan_prn": str(TARGET_PRN),
        "Tracking_1C.state_dump_filename": str(state),
        "Observables.dump": "false",
        "PVT.dump": "false",
    }
    for key, value in values.items():
        prefix = key + "="
        for index, line in enumerate(lines):
            if line.startswith(prefix):
                lines[index] = prefix + value
                break
        else:
            lines.append(prefix + value)
    path = OUT / "state_prn28.conf"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def assess(state: Path) -> dict:
    rows = []
    with state.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            if int(row["prn"]) == TARGET_PRN:
                rows.append({
                    "sample": int(row["tracking_sample_counter"]),
                    "cn0": float(row["cn0_db_hz"]),
                    "lock": float(row["carrier_lock_test"]),
                })
    if not rows:
        raise RuntimeError("no PRN28 state rows")
    first = rows[0]["sample"]
    after = [row for row in rows if row["sample"] >= first + 25_000_000]
    if len(after) < 10_000:
        raise RuntimeError(f"insufficient post-acquisition rows: {len(after)}")
    cn0 = sorted(row["cn0"] for row in after)
    bad = [row["lock"] < 0.5 or row["cn0"] < 30.0 for row in after]
    maximum = current = 0
    for value in bad:
        current = current + 1 if value else 0
        maximum = max(maximum, current)
    locked = 1.0 - sum(bad) / len(after)
    median = cn0[len(cn0) // 2]
    p05 = cn0[max(0, int(0.05 * (len(cn0) - 1)))]
    criteria = {
        "locked_fraction_ge_0p95": locked >= 0.95,
        "cn0_median_ge_35": median >= 35.0,
        "cn0_p05_ge_30": p05 >= 30.0,
        "max_bad_run_ms_le_100": maximum <= 100,
        "post_rows_ge_10000": len(after) >= 10_000,
    }
    return {
        "total_rows": len(rows),
        "post_acquisition_rows": len(after),
        "locked_fraction": locked,
        "cn0_median_db_hz": median,
        "cn0_p05_db_hz": p05,
        "max_consecutive_bad_rows": maximum,
        "criteria": criteria,
        "passed": all(criteria.values()),
    }


def main() -> None:
    if OUT.exists():
        raise RuntimeError(f"refusing to overwrite {OUT}")
    for path in (PROTOCOL, CONFIG, RINEX, SIMULATOR, RECEIVER,
                 STATE_TEMPLATE, SOURCE_GEOMETRY):
        if not path.is_file():
            raise FileNotFoundError(path)
    OUT.mkdir(parents=True)
    TMP.mkdir(parents=True, exist_ok=True)
    if any(TMP.iterdir()):
        raise RuntimeError(f"temporary directory is not empty: {TMP}")
    geometry = extract_geometry()
    iq = TMP / "precheck_iq16.bin"
    state = OUT / "tracking_state_prn28.csv"
    simulator_log = OUT / "simulator.log"
    receiver_log = OUT / "receiver_console.log"
    config = prepare_state(iq, state)
    run([str(SIMULATOR), "--config", str(CONFIG), "--output", str(iq),
         "--noise-seed", str(SEED), "--multi-thread"], simulator_log)
    if iq.stat().st_size != EXPECTED_IQ_BYTES:
        raise RuntimeError(f"unexpected IQ size: {iq.stat().st_size}")
    if "Clipped samples: 0" not in simulator_log.read_text(errors="replace"):
        raise RuntimeError("simulator reported clipping")
    run([str(RECEIVER), f"--config_file={config}", "--logtostderr=1",
         "--colorlogtostderr=1"], receiver_log, OUT)
    metrics = assess(state)
    payload = {
        "status": "prn28_short_tracking_precheck_complete",
        "passed": metrics["passed"], "truth_accessed": False,
        "target_prn": TARGET_PRN, "seed": SEED,
        "protocol_sha256": sha256(PROTOCOL),
        "receiver": str(RECEIVER), "receiver_sha256": sha256(RECEIVER),
        "simulator": str(SIMULATOR), "simulator_sha256": sha256(SIMULATOR),
        "config_sha256": sha256(CONFIG), "state_sha256": sha256(state),
        "geometry": geometry, "quality": metrics,
        "raw_logs": {"simulator": str(simulator_log),
                     "receiver": str(receiver_log), "state": str(state)},
    }
    (OUT / "PRECHECK_REPORT.json").write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    iq.unlink()
    if TMP.exists() and not any(TMP.iterdir()):
        TMP.rmdir()
    print(json.dumps({"passed": payload["passed"], "quality": metrics}, indent=2))
    print(f"REPORT={OUT / 'PRECHECK_REPORT.json'}")


if __name__ == "__main__":
    main()
