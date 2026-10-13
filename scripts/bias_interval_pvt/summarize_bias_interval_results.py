#!/usr/bin/env python3
"""Finalize summaries from already-retained bias interval surfaces and CSVs."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import numpy as np


ROOT = Path("/home/zoxoai/Documents/ChatGPT/paper/PAPER_NEW/results/bias_interval_pvt_20261013")
SURFACES = ROOT / "residual_surfaces"
CHIP_METRES = 299_792_458.0 / 1_023_000.0


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def bool_value(value: str) -> bool:
    return value.strip().lower() in {"true", "1", "yes"}


def uncensored_interval(row: dict[str, str]) -> dict:
    """Recover the profile interval before gate-zeroing, from retained RSS grid."""
    path = Path(row["surface_file"])
    with np.load(path, allow_pickle=False) as data:
        rss = np.asarray(data["rss"], dtype=float)
        tau = np.asarray(data["candidate_tau_chip"], dtype=float)
    minimum = float(np.min(rss))
    kappa = float(row["kappa"])
    near = rss <= minimum * (1.0 + kappa) + 1e-14 * max(minimum, 1.0)
    selected_tau = tau[np.any(near, axis=0)]
    if not selected_tau.size:
        raise RuntimeError(f"empty profile set for {path}")
    tau_zero = float(row["tau_zero_chip"])
    corrections = -(selected_tau - tau_zero) * CHIP_METRES
    return {
        "profile_b_min_m": float(np.min(corrections)),
        "profile_b_max_m": float(np.max(corrections)),
        "profile_width_m": float(np.max(corrections) - np.min(corrections)),
        "profile_tau_min_chip": float(np.min(selected_tau)),
        "profile_tau_max_chip": float(np.max(selected_tau)),
        "profile_candidate_count": int(np.count_nonzero(near)),
    }


def main() -> None:
    interval_rows = read_csv(ROOT / "per_window_bias_intervals.csv")
    coverage_rows = read_csv(ROOT / "short_delay_counterfactual_coverage.csv")
    truth_by_key = {
        (row["seed"], row["cell"], row["window"]): float(row["true_correction_m"])
        for row in coverage_rows
    }
    gate_rows = []
    enriched = []
    for row in interval_rows:
        active = bool_value(row["gate_active"])
        scene = row["scene"]
        injected = scene != "clean"
        gate_rows.append({
            "family": row["family"], "seed": row["seed"], "scene": scene,
            "cell": row["cell"], "window": row["window"],
            "injected_multipath": injected, "gate_active": active,
            "gate_false_positive": bool(active and not injected),
            "gate_false_negative": bool((not active) and injected),
        })
        profile = uncensored_interval(row)
        diagnostic = {**row, **profile}
        if row["family"] == "short_delay":
            true_correction = truth_by_key[(row["seed"], row["cell"], row["window"])]
            diagnostic["true_correction_m"] = true_correction
            diagnostic["profile_covered_ignoring_gate"] = bool(
                profile["profile_b_min_m"] - 1e-9 <= true_correction <= profile["profile_b_max_m"] + 1e-9
            )
        else:
            diagnostic["true_correction_m"] = ""
            diagnostic["profile_covered_ignoring_gate"] = ""
        enriched.append(diagnostic)

    # Preserve the original window-level decisions and add the pre-gate profile
    # interval as a diagnostic; operational b_min/b_max remain [0,0] on gate-off.
    write_csv(ROOT / "gate_window_records.csv", gate_rows)
    write_csv(ROOT / "per_window_bias_intervals_with_profile_diagnostic.csv", enriched)
    write_csv(
        ROOT / "short_delay_ungated_profile_coverage.csv",
        [row for row in enriched if row["family"] == "short_delay"],
    )

    by_cell = {}
    for cell in sorted({row["cell"] for row in coverage_rows}):
        selected = [row for row in coverage_rows if row["cell"] == cell]
        enriched_selected = [row for row in enriched if row["cell"] == cell and row["family"] == "short_delay"]
        by_cell[cell] = {
            "windows": len(selected),
            "coverage_rate": float(np.mean([bool_value(row["covered"]) for row in selected])),
            "median_operational_width_m": float(np.median([float(row["width_m"]) for row in selected])),
            "median_ungated_profile_width_m": float(np.median([float(row["profile_width_m"]) for row in enriched_selected])),
            "gate_active_rate": float(np.mean([bool_value(row["gate_active"]) for row in selected])),
            "ungated_profile_coverage_rate": float(np.mean([
                bool(row["profile_covered_ignoring_gate"]) for row in enriched_selected
            ])),
            "median_true_correction_abs_m": float(np.median([abs(float(row["true_correction_m"])) for row in selected])),
        }

    by_delay = {}
    for delay, prefix in ((0.05, "d005"), (0.10, "d010")):
        selected = [row for row in coverage_rows if row["cell"].startswith(prefix)]
        enriched_selected = [row for row in enriched if row["family"] == "short_delay" and row["cell"].startswith(prefix)]
        by_delay[str(delay)] = {
            "windows": len(selected),
            "coverage_rate": float(np.mean([bool_value(row["covered"]) for row in selected])),
            "median_operational_width_m": float(np.median([float(row["width_m"]) for row in selected])),
            "median_ungated_profile_width_m": float(np.median([float(row["profile_width_m"]) for row in enriched_selected])),
            "gate_active_rate": float(np.mean([bool_value(row["gate_active"]) for row in selected])),
            "ungated_profile_coverage_rate": float(np.mean([
                bool(row["profile_covered_ignoring_gate"]) for row in enriched_selected
            ])),
        }

    gate_summary = {}
    for family, scenes in (("long_0p8", ("clean", "amp_0p05", "amp_0p10")), ("short_delay", ("injected",))):
        for scene in scenes:
            selected = [row for row in gate_rows if row["family"] == family and row["scene"] == scene]
            if not selected:
                continue
            active_rate = float(np.mean([bool(row["gate_active"]) for row in selected]))
            injected = bool(selected[0]["injected_multipath"])
            gate_summary[f"{family}/{scene}"] = {
                "windows": len(selected), "active_rate": active_rate,
                "false_positive_rate": active_rate if not injected else None,
                "false_negative_rate": 1.0 - active_rate if injected else None,
            }

    short_coverage = float(np.mean([bool_value(row["covered"]) for row in coverage_rows]))
    kappas = {
        family: sorted({float(row["kappa"]) for row in interval_rows if row["family"] == family})
        for family in ("long_0p8", "short_delay")
    }
    tau_zero = {
        family: sorted({float(row["tau_zero_chip"]) for row in interval_rows if row["family"] == family})
        for family in ("long_0p8", "short_delay")
    }
    protocol = ROOT / "FROZEN_REANALYSIS_RULES_20261013.md"
    summary = {
        "status": "stage1_and_stage2_existing_data_reanalysis_complete",
        "analysis_class": "exploratory_post_hoc_reanalysis",
        "data_generated": False,
        "protocol_sha256": sha256(protocol),
        "grid": {
            "candidate_pairs_tau_separation": 2419,
            "rho_values": 113,
            "residual_surface_files": len(list(SURFACES.rglob("*.npz"))),
            "grid_points_per_window": 2419 * 113,
        },
        "calibration": {"kappa_by_family": kappas, "tau_zero_values_by_family": tau_zero},
        "long_0p8": {
            "windows": sum(row["family"] == "long_0p8" for row in interval_rows),
            "counterfactual_coverage": None,
            "coverage_status": "not scoreable: retained batch lacks same-seed same-second clean RINEX counterfactual",
        },
        "short_delay": {
            "windows": len(coverage_rows),
            "coverage_rate": short_coverage,
            "median_operational_width_m": float(np.median([float(row["width_m"]) for row in coverage_rows])),
            "median_ungated_profile_width_m": float(np.median([
                float(row["profile_width_m"]) for row in enriched if row["family"] == "short_delay"
            ])),
            "ungated_profile_coverage_rate": float(np.mean([
                bool(row["profile_covered_ignoring_gate"])
                for row in enriched if row["family"] == "short_delay"
            ])),
            "coverage_by_cell": by_cell,
            "coverage_by_delay_chip": by_delay,
        },
        "gate_rates": gate_summary,
        "stage2_decision": "PASS_TO_PVT" if short_coverage >= 0.90 else "STOP_BEFORE_PVT",
    }
    (ROOT / "SUMMARY.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    md = [
        "# Bias-interval reanalysis results",
        "",
        "Classification: exploratory/post-hoc; all inputs are retained CSV/RINEX/summary data; no IQ was generated.",
        "",
        f"- Full grid residuals retained: {summary['grid']['residual_surface_files']} windows × {summary['grid']['grid_points_per_window']:,} candidates.",
        f"- κ (0.8 chip): {kappas['long_0p8']}; κ (short delay): {kappas['short_delay']}.",
        f"- Comparable v3 clean gate false-positive rate: {gate_summary['long_0p8/clean']['false_positive_rate']:.3f}.",
        "- Comparable v3 0.8-chip counterfactual coverage: not scoreable; no paired clean RINEX exists for the same seed and second.",
        f"- Short-delay counterfactual coverage: {short_coverage:.3f} (n={len(coverage_rows)} windows); median operational width {summary['short_delay']['median_operational_width_m']:.3f} m.",
        f"- Ignoring the gate, profile-envelope coverage: {summary['short_delay']['ungated_profile_coverage_rate']:.3f}; median profile width {summary['short_delay']['median_ungated_profile_width_m']:.3f} m.",
        f"- Short-delay gate false-negative rate: {gate_summary['short_delay/injected']['false_negative_rate']:.3f}.",
        f"- Stage-2 decision under the frozen ≥90% rule: **{summary['stage2_decision']}**.",
        "",
        "| Cell | Operational coverage | Gate active | Operational width (median m) | Ungated profile coverage | Ungated width (median m) |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for cell, item in by_cell.items():
        md.append(f"| {cell} | {item['coverage_rate']:.3f} | {item['gate_active_rate']:.3f} | {item['median_operational_width_m']:.3f} | {item['ungated_profile_coverage_rate']:.3f} | {item['median_ungated_profile_width_m']:.3f} |")
    md += [
        "",
        "Since step 2 fails, the bounded-PVT and M1–M8 comparisons were not run. The retained surface files allow the interval calculations to be audited without repeating the fits.",
        "",
        "Files: `per_window_bias_intervals.csv`, `per_window_bias_intervals_with_profile_diagnostic.csv`, `short_delay_counterfactual_coverage.csv`, `short_delay_ungated_profile_coverage.csv`, and `residual_surfaces/`.",
    ]
    (ROOT / "REPORT.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": summary["status"], "short_coverage": short_coverage,
        "short_gate_active_rate": gate_summary["short_delay/injected"]["active_rate"],
        "median_short_ungated_width_m": summary["short_delay"]["median_ungated_profile_width_m"],
        "stage2_decision": summary["stage2_decision"],
    }, indent=2))


if __name__ == "__main__":
    main()
