"""Frozen exploratory replay: one retained seed, five scenes, no IQ generation."""
import csv
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import time

import numpy as np
from scripts.baselines.correlator_baselines import (
    BaselineConfig, EpochScan, eml_estimate, double_delta_estimate,
)
from scripts.gss_multipath.gss import Config, SampledCorrelation, estimate

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/gss_replay_20261020"
SEED = ROOT / "results/prn09_phase_amp_delay_grid_20261019/seed_2026101922"
SCENES = ["clean", "d0p05_a0p20_p000", "d0p05_a0p20_p090",
          "d0p05_a0p20_p180", "d0p05_a0p50_p090"]
WINDOWS = list(range(27, 80))
METRES = 299792458.0 / 1023000


def process(path, destination):
    estimates = []
    config = BaselineConfig()
    with path.open() as stream:
        groups = itertools.groupby(csv.DictReader(stream),
                                   lambda r: (int(r["window"]), int(r["epoch"])))
        for (window, epoch), grouped in groups:
            rows = list(grouped)
            x = -np.asarray([float(r["delay_chip"]) for r in rows])[::-1]
            output = {"window": window, "epoch": epoch,
                      "time_ms": float(rows[0]["time_ms"])}
            for arm in ("primary", "mixed"):
                z = np.asarray([complex(float(r[f"{arm}_i"]), float(r[f"{arm}_q"]))
                                for r in rows])[::-1]
                correlation = SampledCorrelation(x, z)
                correlation.validate_support(Config())
                result = estimate(correlation)
                scan = EpochScan(window, epoch, window, x, z, {})
                output[f"{arm}_gss"] = result["estimate_chip"]
                output[f"{arm}_eml"] = eml_estimate(scan, config)["tau_chip"]
                output[f"{arm}_dd"] = double_delta_estimate(scan, config)["tau_chip"]
                output[f"{arm}_gss_alias_count"] = sum(
                    r["one_chip_alias_rejected"] for r in result["rounds"])
                output[f"{arm}_gss_round_spread"] = float(np.ptp(
                    [r["candidate_chip"] for r in result["rounds"]]))
            estimates.append(output)
    with destination.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(estimates[0]))
        writer.writeheader()
        writer.writerows(estimates)
    return estimates


def main():
    OUT.mkdir(exist_ok=True, parents=True)
    protocol = {
        "type": "exploratory_offline_interpolated_GSS_not_closed_loop",
        "seed": 2026101922, "scenes": SCENES, "windows": WINDOWS,
        "calibration_windows": list(range(27, 40)),
        "evaluation_windows": list(range(40, 80)),
        "coherent_epochs_per_window": 100, "window_reducer": "arithmetic_mean",
        "reference": "unit_triangle", "gate": "none; standalone paper estimator",
        "comparators": "existing EML 0.1 / DD inner 0.1 outer 0.5",
        "coordinate": "physical_delay = -CUDA_receiver_tap_shift",
        "scanner_limit": "same original replay; front-end FIR not replayed",
        "scope": "single retained seed, posthoc data, no confirmatory inference",
    }
    frozen = OUT / "REPLAY_PROTOCOL.json"
    if frozen.exists():
        if json.loads(frozen.read_text()) != protocol:
            raise RuntimeError("protocol changed")
    else:
        frozen.write_text(json.dumps(protocol, indent=2))
    scanner = OUT / "build/gss_wide_scan"
    for scene in SCENES:
        all_rows = []
        for start in range(0, len(WINDOWS), 13):
            chunk = WINDOWS[start:start+13]
            directory = OUT / scene / f"chunk_{start//13:02d}"
            directory.mkdir(parents=True, exist_ok=True)
            epochs = directory / "epochs.csv"
            result = directory / "estimates.csv"
            if result.exists():
                with result.open() as stream:
                    all_rows.extend([{k: float(v) for k, v in r.items()}
                                     for r in csv.DictReader(stream)])
                continue
            command = [
                str(scanner), "--prn", "9", "--state", str(SEED/scene/"tracking_state.csv"),
                "--primary", str(SEED/"clean/signal.iq16"),
                "--secondary", str(SEED/"clean/signal.iq16"),
                "--mixed", str(SEED/scene/"signal.iq16"),
                "--output-dir", str(directory), "--epoch-output", str(epochs),
                "--epochs-per-window", "100", "--max-windows", str(len(chunk)),
                "--grid-ms", "1000", "--skip-ms", str(chunk[0]*1000),
                "--fs-hz", "25000000", "--if-hz", "-4580000",
            ]
            started = time.monotonic()
            print("SCAN", scene, chunk[0], chunk[-1], flush=True)
            if not (directory/"SCAN_COMPLETE.json").exists():
                with (directory/"scanner.log").open("w") as log:
                    subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
                (directory/"SCAN_COMPLETE.json").write_text(json.dumps({
                    "command": command, "elapsed_seconds": time.monotonic()-started,
                    "binary_sha256": hashlib.sha256(scanner.read_bytes()).hexdigest(),
                }, indent=2))
            all_rows.extend(process(epochs, result))
            print("DONE", scene, chunk[0], round(time.monotonic()-started, 1), flush=True)
        if len(all_rows) != len(WINDOWS)*100:
            raise RuntimeError("incomplete epoch count")
        summary = []
        for window in WINDOWS:
            rows = [r for r in all_rows if int(r["window"]) == window]
            entry = {"window": window}
            for arm in ("primary", "mixed"):
                for method in ("gss", "eml", "dd"):
                    entry[f"{arm}_{method}"] = float(np.mean(
                        [r[f"{arm}_{method}"] for r in rows]))
            summary.append(entry)
        (OUT/scene/"window_estimates.json").write_text(json.dumps(summary, indent=2))
    report = {}
    for scene in SCENES:
        rows = json.loads((OUT/scene/"window_estimates.json").read_text())
        rows = [r for r in rows if r["window"] >= 40]
        report[scene] = {}
        for method in ("gss", "eml", "dd"):
            difference = np.asarray([r[f"mixed_{method}"]-r[f"primary_{method}"]
                                     for r in rows])*METRES
            report[scene][method] = {
                "paired_clean_displacement_rmse_m": float(np.sqrt(np.mean(difference**2))),
                "paired_clean_displacement_mean_m": float(np.mean(difference)),
                "windows": len(rows),
            }
    (OUT/"CORRELATOR_DIAGNOSTIC.json").write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
