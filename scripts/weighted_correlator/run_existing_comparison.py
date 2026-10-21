"""Compare a weighted early/late correction with existing PRN9 observations.

This is a posthoc, one-seed common-observation pilot. It does not tune weights
on the evaluation windows and does not alter the GNSS-SDR receiver binary.
"""
import csv
import hashlib
import itertools
import json
from datetime import timedelta
from pathlib import Path

import numpy as np

from scripts.gss_multipath.replay import OUT as REPLAY_OUT, SEED, SCENES, METRES
from scripts.weighted_correlator.core import Configuration, estimate
from scripts.prn09_phase_delay_amplitude_grid import complete_seed22_gated_pvt as helper


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/weighted_correlator_comparison_20261021"
CAL_WINDOWS = list(range(27, 40))
EVAL_WINDOWS = list(range(40, 80))
CONFIG = Configuration(spacings=(0.10, 0.20), weights=(1.0, -0.5),
                        radius=0.15, step=0.001)


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def estimate_file(path):
    outputs = []
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        groups = itertools.groupby(reader, key=lambda r: (int(r["window"]), int(r["epoch"])))
        for (window, epoch), group in groups:
            taps = list(group)
            if len(taps) != 401:
                raise ValueError(f"{path}: window {window} epoch {epoch} has {len(taps)} taps")
            taps.sort(key=lambda r: float(r["delay_chip"]))
            x = -np.asarray([float(r["delay_chip"]) for r in taps])[::-1]
            primary = np.asarray([complex(float(r["primary_i"]), float(r["primary_q"]))
                                  for r in taps])[::-1]
            mixed = np.asarray([complex(float(r["mixed_i"]), float(r["mixed_q"]))
                                for r in taps])[::-1]
            outputs.append({"window": window, "epoch": epoch,
                            "primary": estimate(x, primary, reference=0.0, config=CONFIG),
                            "mixed": estimate(x, mixed, reference=0.0, config=CONFIG)})
    return outputs


def main():
    if OUT.exists():
        raise FileExistsError(f"refusing to overwrite {OUT}")
    OUT.mkdir(parents=True)
    epoch_rows = []
    per_scene = {}
    input_hashes = {}
    for scene in SCENES:
        chunks = sorted((REPLAY_OUT / scene).glob("chunk_*/epochs.csv"))
        if not chunks:
            raise FileNotFoundError(f"no retained correlator scans for {scene}")
        all_rows = []
        for path in chunks:
            input_hashes[str(path)] = digest(path)
            all_rows.extend(estimate_file(path))
        grouped = {}
        for row in all_rows:
            grouped.setdefault(row["window"], []).append(row)
            for arm in ("primary", "mixed"):
                result = row[arm]
                epoch_rows.append({"scene": scene, "window": row["window"],
                                   "epoch": row["epoch"], "arm": arm,
                                   "valid": result["valid"],
                                   "tau_chip": result["estimate_chip"],
                                   "reason": result["reason"]})
        summaries = {}
        for window, rows in grouped.items():
            record = {"window": window, "epoch_count": len(rows)}
            for arm in ("primary", "mixed"):
                valid = [r[arm]["estimate_chip"] for r in rows if r[arm]["valid"]]
                record[f"{arm}_valid_count"] = len(valid)
                record[f"{arm}_tau_chip"] = float(np.mean(valid)) if valid else None
            summaries[window] = record
        per_scene[scene] = summaries
    with (OUT / "EPOCH_ESTIMATES.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(epoch_rows[0]))
        writer.writeheader(); writer.writerows(epoch_rows)

    clean = per_scene["clean"]
    calibration_values = [clean[w]["mixed_tau_chip"] for w in CAL_WINDOWS
                          if clean[w]["mixed_valid_count"] == clean[w]["epoch_count"]]
    if len(calibration_values) != len(CAL_WINDOWS):
        raise RuntimeError("weighted discriminator did not return 100/100 valid epochs in every clean calibration window")
    zero = float(np.median(calibration_values))
    calibration = {"method": "median of clean mixed weighted-discriminator window means",
                   "windows": CAL_WINDOWS, "zero_chip": zero,
                   "zero_m": zero * METRES, "weights": CONFIG.weights,
                   "spacings_chip": CONFIG.spacings}
    (OUT / "CALIBRATION.json").write_text(json.dumps(calibration, indent=2) + "\n")

    summary_rows = []
    for scene, windows in per_scene.items():
        for window, record in windows.items():
            if record["mixed_tau_chip"] is None:
                correction = 0.0
            elif scene == "clean" and window in CAL_WINDOWS:
                correction = 0.0
            elif record["mixed_valid_count"] != record["epoch_count"]:
                # Fail closed: do not apply a correction if the estimator was
                # invalid for any epoch in this one-second tracking window.
                correction = 0.0
            else:
                correction = -(record["mixed_tau_chip"] - zero) * METRES
            summary_rows.append({"scene": scene, **record,
                                 "weighted_correction_m": correction})
    with (OUT / "WINDOW_ESTIMATES.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(summary_rows[0]))
        writer.writeheader(); writer.writerows(summary_rows)
    for scene, windows in per_scene.items():
        (OUT / f"{scene}_window_estimates.json").write_text(
            json.dumps([windows[w] for w in sorted(windows)], indent=2) + "\n")

    protocol = {"type": "posthoc_existing_data_single_seed",
                "source": "retained PRN9 seed2026101922 401-tap per-epoch scans",
                "calibration_windows": CAL_WINDOWS, "evaluation_windows": EVAL_WINDOWS,
                "discriminator": "coherent weighted early-late sum, spacings .10/.20 chip, weights 1/-.5",
                "method_scope": "weighted correlator prototype; no adaptive MAA lookup; no tuned weights",
                "zero_reference": calibration,
                "invalid_window_policy": "zero correction unless every 1ms epoch has a valid root",
                "inputs_sha256": input_hashes,
                "truth_used_for_estimation": False,
                "caveat": "offline observation correction, not closed-loop tracking"
                }
    (OUT / "PROTOCOL_AND_INPUTS.json").write_text(json.dumps(protocol, indent=2) + "\n")

    # Create weighted-correction observation branches and solve with exactly
    # the receiver/PVT helper used by the retained common-observation results.
    helper.SCENES = {scene: SEED / scene for scene in SCENES}
    pvt = {}
    for scene in SCENES:
        directory = OUT / scene / "pvt_weighted"
        directory.mkdir(parents=True)
        schedule = directory / "corrections.csv"
        with schedule.open("w", newline="", encoding="utf-8") as f:
            fields = ("start_sample", "end_sample", "prn", "correction_m")
            writer = csv.DictWriter(f, fieldnames=fields); writer.writeheader()
            for window in range(80):
                estimate_row = per_scene[scene].get(window)
                correction = 0.0
                if estimate_row is not None and window >= 40:
                    if estimate_row["mixed_valid_count"] == estimate_row["epoch_count"]:
                        correction = -(estimate_row["mixed_tau_chip"] - zero) * METRES
                writer.writerow({"start_sample": window * 25_000_000,
                                 "end_sample": (window+1)*25_000_000-1,
                                 "prn": 9, "correction_m": correction})
        paths = helper.build_rinex_branches(
            scene, directory, schedule,
            helper.SIGNAL_START_UTC + timedelta(seconds=40),
            helper.SIGNAL_START_UTC + timedelta(seconds=80))
        solutions = helper.solve_pair(scene, paths, directory)
        truth = (37.352721, -121.915773, 20.0)
        parsed = {k: helper.parse_solution(Path(v["path"])) for k,v in solutions.items()}
        expected = paths["score_seconds"]
        available = sorted(set(expected) & set(parsed["gated"]) & set(parsed["single"]))
        errors = {}
        for branch in ("single", "gated"):
            enu = np.asarray([helper.enu_error(parsed[branch][s]["latitude_deg"],
                                parsed[branch][s]["longitude_deg"],
                                parsed[branch][s]["height_m"], *truth) for s in available])
            d3 = np.linalg.norm(enu, axis=1)
            errors[branch] = {"3d_rmse_m": float(np.sqrt(np.mean(d3**2))),
                              "horizontal_rmse_m": float(np.sqrt(np.mean(np.sum(enu[:,:2]**2,axis=1)))),
                              "vertical_rmse_m": float(np.sqrt(np.mean(enu[:,2]**2))),
                              "p95_m": float(np.percentile(d3,95)),
                              "gross_gt10m": float(np.mean(d3>10)),
                              "common_epochs": len(available),
                              "missing_epochs": sorted(set(expected)-set(available))}
        pvt[scene] = {"paths": paths, "solutions": solutions, "scores": errors}
    (OUT / "WEIGHTED_PVT.json").write_text(json.dumps(pvt, indent=2) + "\n")

    # Join previously computed original/EML/DD branches from the exact same
    # retained-observation comparison and report the new weighted branch.
    previous = json.loads((REPLAY_OUT / "PVT_SCORING.json").read_text())
    rows = []
    for scene in SCENES:
        refs = previous["scenes"][scene]["methods"]
        for method in ("original", "eml", "dd"):
            rows.append({"scene": scene, "method": method,
                         **{k: refs[method][k] for k in
                            ("three_d_rmse_m","horizontal_rmse_m","vertical_rmse_m","p95_m","gross_gt10m","common_epochs")
                            if k in refs[method]}})
        weighted_score = pvt[scene]["scores"]["gated"]
        rows.append({"scene": scene, "method": "weighted_EL",
                     "three_d_rmse_m": weighted_score["3d_rmse_m"],
                     "horizontal_rmse_m": weighted_score["horizontal_rmse_m"],
                     "vertical_rmse_m": weighted_score["vertical_rmse_m"],
                     "p95_m": weighted_score["p95_m"],
                     "gross_gt10m": weighted_score["gross_gt10m"],
                     "common_epochs": weighted_score["common_epochs"]})
    with (OUT / "COMPARISON.csv").open("w", newline="", encoding="utf-8") as f:
        fields = list(dict.fromkeys(key for row in rows for key in row))
        writer = csv.DictWriter(f, fieldnames=fields); writer.writeheader(); writer.writerows(rows)
    print(json.dumps({"output": str(OUT), "zero": calibration,
                      "comparison": rows}, indent=2))


if __name__ == "__main__":
    main()
