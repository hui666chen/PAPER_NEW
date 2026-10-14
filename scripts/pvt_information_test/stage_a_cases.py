#!/usr/bin/env python3
"""Stage A: per-satellite semi-synthetic cases (no IQ generation, no PVT).

Known-truth reflections are added to REAL clean 401-tap observations
(100-epoch means and epoch bootstraps of the PRN09 clean calibration windows),
so the noise, receiver-response mismatch and affine background are real while
the reflection parameters are known.

For every clean/injected observation the script records
  * the Double-Delta (DD) and EML residual pseudorange error (truth),
  * the two-path profile interval for the DD error   [lo, hi]  (proposed),
  * the single-vs-two-path detection statistic T     (for exclusion baselines).

Calibration/test split (frozen): within each cell the 8 clean windows are
split by position, even -> calibration, odd -> test.  The DD/EML zero points,
the detection threshold and the interval tolerance kappa are fixed on the
calibration windows only, then every test row is scored once.

Usage on the experiment PC (paths default to the existing layout):
    python3 scripts/pvt_information_test/stage_a_cases.py --workers 4
Self-test with a synthetic response (no data needed):
    python3 scripts/pvt_information_test/stage_a_cases.py --self-test
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path
import sys
import time

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REPO / "scripts" / "baselines"))

from common import CHIP_METRES, StableProfile, quantile_higher, read_epoch_matrix, write_csv  # noqa: E402
from correlator_baselines import (  # noqa: E402
    BaselineConfig, EpochScan, double_delta_estimate, eml_estimate)

DEFAULT_PAPER = REPO.parent
BANKS = {"d005_a025": 39, "d005_a050": 53, "d010_a025": 67, "d010_a050": 81}
CLEAN_WINDOWS = {
    "d005_a025": tuple(range(36, 44)),
    "d005_a050": tuple(range(50, 58)),
    "d010_a025": tuple(range(64, 72)),
    "d010_a050": tuple(range(78, 86)),
}
# Frozen injection design (operating envelope for both calibration and test).
DELAYS_CHIP = (0.02, 0.05, 0.10, 0.20)
AMPLITUDES = (0.10, 0.25, 0.50)
PHASES_DEG = tuple(range(0, 360, 45))
KAPPA_QUANTILE = 0.95
DETECTION_QUANTILE = 0.99
SEED = 2026101420


# ---------------------------------------------------------------- data access
def load_real_cell(cell: str, data_dir: Path, paper: Path):
    sys.path.insert(0, str(paper))
    sys.path.insert(0, str(paper / "PAPER_NEW" / "scripts"))
    sys.path.insert(0, str(paper / "gps_l1_prn09_shared_ratio_v2"))
    sys.path.insert(0, str(REPO / "scripts"))
    try:
        from receiver_chain_response_control import ResponseSpec, load_response_model
    except ImportError as exc:  # pragma: no cover - depends on the local PC
        raise SystemExit(
            f"Cannot import receiver_chain_response_control ({exc}). Run on the experiment PC "
            f"where scripts/short_delay_detectability/refine_bias_intervals.py runs, or pass --paper.")
    spec = ResponseSpec(prn=9)
    model = load_response_model(data_dir / "response_banks" / f"bank_{BANKS[cell]}.npz",
                                expected_spec=spec)
    epochs = read_epoch_matrix(data_dir / "open_loop" / cell / "prn09_epoch_responses.csv",
                               CLEAN_WINDOWS[cell])
    return model, spec.tap_grid.astype(np.float64), epochs


class SyntheticModel:
    """Band-limited triangle used only by --self-test."""

    def __init__(self, tap: np.ndarray):
        self.tap = tap
        self.kernel = np.linspace(-0.08, 0.08, 9)

    def exact(self, delay_chip):
        d = np.atleast_1d(np.asarray(delay_chip, dtype=float))
        x = self.tap[None, :, None] + d[:, None, None] + self.kernel[None, None, :]
        return (np.maximum(0.0, 1.0 - np.abs(x)).mean(axis=2) * np.exp(0.3j)).astype(np.complex128)


def load_synthetic_cell(cell: str, rng: np.random.Generator):
    tap = -2.0 + 0.01 * np.arange(401)
    model = SyntheticModel(tap)
    epochs = {}
    for w in CLEAN_WINDOWS[cell]:
        tau_w = rng.normal(0.0, 0.004)
        clean = 1000.0 * model.exact([-tau_w])[0]
        affine = (3.0 + 2.0j) + (1.5 - 1.0j) * tap
        noise = (rng.normal(size=(100, 401)) + 1j * rng.normal(size=(100, 401))) * 400.0
        epochs[w] = clean[None, :] + affine[None, :] + noise
    return model, tap, epochs


# ---------------------------------------------------------------- measurement
class Measurer:
    def __init__(self, model, tap, step):
        self.model = model
        self.tap = tap
        self.profile = StableProfile(model, tap, step=step)
        self.cfg = BaselineConfig()

    def __call__(self, y: np.ndarray) -> dict:
        prof = self.profile.profile(y)
        tau1, gain = self.profile.onepath_fit(prof)
        scan = EpochScan(0, 0, 0, self.tap, y, {})
        return {
            "prof": prof,
            "tau1": tau1,
            "gain": gain,
            "dd": float(double_delta_estimate(scan, self.cfg)["tau_chip"]),
            "eml": float(eml_estimate(scan, self.cfg)["tau_chip"]),
            "T": StableProfile.detection_stat(prof),
        }

    def reflection(self, tau1: float, gain: complex, delay: float, amp: float, phase_deg: float):
        rho = amp * np.exp(1j * np.deg2rad(phase_deg))
        return gain * rho * np.asarray(self.model.exact(np.asarray([delay - tau1])))[0]


def bootstrap_means(epochs: np.ndarray, count: int, rng: np.random.Generator) -> list[np.ndarray]:
    means = [epochs.mean(axis=0)]
    for _ in range(count):
        means.append(epochs[rng.integers(0, len(epochs), size=len(epochs))].mean(axis=0))
    return means


def designs():
    for delay in DELAYS_CHIP:
        for amp in AMPLITUDES:
            for phase in PHASES_DEG:
                yield delay, amp, phase


def process_cell(job: dict) -> dict:
    t0 = time.time()
    cell = job["cell"]
    rng = np.random.default_rng(job["seed"])
    if job["self_test"]:
        model, tap, epochs = load_synthetic_cell(cell, rng)
    else:
        model, tap, epochs = load_real_cell(cell, Path(job["data_dir"]), Path(job["paper"]))
    meas = Measurer(model, tap, job["step"])
    prof_obj = meas.profile
    windows = sorted(epochs)
    cal_w, test_w = windows[0::2], windows[1::2]

    # 1) Clean calibration: zero points and detection threshold.
    cal_clean = []
    for w in cal_w:
        for b, y in enumerate(bootstrap_means(epochs[w], job["cal_boot"], rng)):
            m = meas(y)
            m.update(window=w, boot=b, y=y)
            cal_clean.append(m)
    c_dd = float(np.median([m["dd"] - m["tau1"] for m in cal_clean]))
    c_eml = float(np.median([m["eml"] - m["tau1"] for m in cal_clean]))
    t_thr = quantile_higher([m["T"] for m in cal_clean], DETECTION_QUANTILE)

    # 2) Injected calibration: kappa so that the retained hull contains the
    #    paired clean direct delay in KAPPA_QUANTILE of calibration injections.
    required = []
    for m in cal_clean:
        if m["boot"] > job["cal_inject_boot"]:
            continue
        required.append(prof_obj.required_kappa(m["prof"], m["tau1"]))
        for delay, amp, phase in designs():
            ymix = m["y"] + meas.reflection(m["tau1"], m["gain"], delay, amp, phase)
            prof = prof_obj.profile(ymix)
            required.append(prof_obj.required_kappa(prof, m["tau1"]))
    required = np.asarray(required)
    kappa = quantile_higher(required[np.isfinite(required)], KAPPA_QUANTILE)

    # 3) Test windows: score every clean and injected observation once.
    half_m = 0.5 * prof_obj.step * CHIP_METRES
    edge = (float(prof_obj.tau[0]) + prof_obj.step, float(prof_obj.tau[-1]) - prof_obj.step)
    rows = []

    def row(base, mix, w, b, delay, amp, phase):
        s_min, s_max, tau_hat = prof_obj.interval(mix["prof"], kappa)
        dd = mix["dd"]
        e_dd = (base["tau1"] - dd + c_dd) * CHIP_METRES
        lo = (s_min - dd + c_dd) * CHIP_METRES
        hi = (s_max - dd + c_dd) * CHIP_METRES
        return {
            "cell": cell, "window": w, "boot": b,
            "scene": "clean" if amp == 0 else "injected",
            "delay_chip": delay, "amplitude": amp, "phase_deg": phase,
            "tau1_clean_chip": base["tau1"], "tau1_mix_chip": mix["tau1"],
            "tau_dd_chip": dd, "tau_eml_chip": mix["eml"],
            "tau_profile_hat_chip": tau_hat, "s_min_chip": s_min, "s_max_chip": s_max,
            "e_dd_m": e_dd,
            "e_eml_m": (base["tau1"] - mix["eml"] + c_eml) * CHIP_METRES,
            "e_profile_point_m": (base["tau1"] - tau_hat) * CHIP_METRES,
            "e_onepath_m": (base["tau1"] - mix["tau1"]) * CHIP_METRES,
            "dd_interval_lo_m": lo, "dd_interval_hi_m": hi, "dd_interval_width_m": hi - lo,
            "covered": bool(lo - half_m <= e_dd <= hi + half_m),
            "interval_at_grid_edge": bool(s_min <= edge[0] or s_max >= edge[1]),
            "T": mix["T"], "flagged": bool(mix["T"] > t_thr),
        }

    for w in test_w:
        for b, y in enumerate(bootstrap_means(epochs[w], job["test_boot"], rng)):
            base = meas(y)
            rows.append(row(base, base, w, b, 0.0, 0.0, 0.0))
            for delay, amp, phase in designs():
                mix = meas(y + meas.reflection(base["tau1"], base["gain"], delay, amp, phase))
                rows.append(row(base, mix, w, b, delay, amp, phase))

    calib = {
        "cell": cell, "calibration_windows": cal_w, "test_windows": test_w,
        "c_dd_chip": c_dd, "c_eml_chip": c_eml, "detection_threshold_T": t_thr,
        "kappa": kappa, "kappa_quantile": KAPPA_QUANTILE,
        "kappa_calibration_count": int(np.isfinite(required).sum()),
        "clean_calibration_observations": len(cal_clean),
        "profile_step_chip": prof_obj.step, "seconds": round(time.time() - t0, 1),
    }
    print(f"[stage A] {cell}: kappa={kappa:.4g} T_thr={t_thr:.4g} rows={len(rows)} "
          f"({calib['seconds']} s)", flush=True)
    return {"rows": rows, "calib": calib}


# ---------------------------------------------------------------- summary
def summarize(rows: list[dict]) -> dict:
    out = {"cells": {}}
    inj = [r for r in rows if r["scene"] == "injected"]
    cln = [r for r in rows if r["scene"] == "clean"]
    out["clean_coverage"] = float(np.mean([r["covered"] for r in cln]))
    out["injected_coverage"] = float(np.mean([r["covered"] for r in inj]))
    out["clean_flag_rate"] = float(np.mean([r["flagged"] for r in cln]))
    table = []
    groups = {}
    for r in inj:
        groups.setdefault((r["delay_chip"], r["amplitude"]), []).append(r)
    for (delay, amp), g in sorted(groups.items()):
        lo = np.asarray([r["dd_interval_lo_m"] for r in g], dtype=float)
        hi = np.asarray([r["dd_interval_hi_m"] for r in g], dtype=float)
        e = np.asarray([r["e_dd_m"] for r in g], dtype=float)
        big = np.abs(e) > 1.0
        mid = 0.5 * (lo + hi)
        table.append({
            "delay_chip": delay, "amplitude": amp, "n": len(g),
            "coverage": float(np.mean([r["covered"] for r in g])),
            "median_width_m": float(np.median([r["dd_interval_width_m"] for r in g])),
            "one_sided_fraction": float(np.mean((lo > 0) | (hi < 0))),
            "sign_ok_when_abs_err_gt_1m": float(np.mean(np.sign(mid[big]) == np.sign(e[big]))) if big.any() else float("nan"),
            "rms_dd_m": float(np.sqrt((e ** 2).mean())),
            "rms_eml_m": float(np.sqrt(np.mean(np.asarray([r["e_eml_m"] for r in g], dtype=float) ** 2))),
            "rms_profile_point_m": float(np.sqrt(np.mean(np.asarray([r["e_profile_point_m"] for r in g], dtype=float) ** 2))),
            "rms_after_midpoint_m": float(np.sqrt(np.mean((e - mid) ** 2))),
            "detect_rate": float(np.mean([r["flagged"] for r in g])),
            "grid_edge_rate": float(np.mean([r["interval_at_grid_edge"] for r in g])),
        })
    out["by_delay_amplitude"] = table
    big_cells = [r for r in table if r["delay_chip"] >= 0.05]
    worst = min(r["coverage"] for r in big_cells) if big_cells else float("nan")
    ok = (out["injected_coverage"] >= 0.90 and out["clean_coverage"] >= 0.90 and worst >= 0.80)
    out["worst_cell_coverage_delay_ge_0p05"] = worst
    out["stage_a_decision"] = "STAGE_A_PASS" if ok else "STAGE_A_STOP_INTERVALS_UNSAFE"
    out["stage_a_rule"] = ("injected coverage >= 0.90, clean coverage >= 0.90, "
                           "and every (delay>=0.05, amplitude) cell >= 0.80")
    return out


def write_summary_md(path: Path, summary: dict, calib: list[dict]) -> None:
    lines = ["# Stage A: per-satellite semi-synthetic interval check", "",
             f"Decision: **{summary['stage_a_decision']}**  ",
             f"Rule: {summary['stage_a_rule']}", "",
             f"- Clean coverage: {summary['clean_coverage']:.3f}; clean flag rate: {summary['clean_flag_rate']:.3f}",
             f"- Injected coverage: {summary['injected_coverage']:.3f}; worst cell (delay>=0.05): "
             f"{summary['worst_cell_coverage_delay_ge_0p05']:.3f}", "",
             "| delay (chip) | amp | n | coverage | median width (m) | one-sided | sign ok (|e|>1 m) | "
             "RMS DD (m) | RMS EML (m) | RMS profile point (m) | RMS after midpoint (m) | detect | grid edge |",
             "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for r in summary["by_delay_amplitude"]:
        lines.append(
            f"| {r['delay_chip']:.2f} | {r['amplitude']:.2f} | {r['n']} | {r['coverage']:.3f} | "
            f"{r['median_width_m']:.2f} | {r['one_sided_fraction']:.3f} | {r['sign_ok_when_abs_err_gt_1m']:.3f} | "
            f"{r['rms_dd_m']:.2f} | {r['rms_eml_m']:.2f} | {r['rms_profile_point_m']:.2f} | "
            f"{r['rms_after_midpoint_m']:.2f} | {r['detect_rate']:.3f} | {r['grid_edge_rate']:.3f} |")
    lines += ["", "## Calibration (calibration windows only)", "",
              "| cell | kappa | T threshold | c_DD (chip) | c_EML (chip) |", "|---|---:|---:|---:|---:|"]
    for c in calib:
        lines.append(f"| {c['cell']} | {c['kappa']:.4g} | {c['detection_threshold_T']:.4g} | "
                     f"{c['c_dd_chip']:.5f} | {c['c_eml_chip']:.5f} |")
    lines += ["", "Noise realisations within a window are epoch bootstraps of the same 100 epochs; "
              "the independent unit is the clean window (4 test windows per cell). "
              "Truth for the direct-path delay is the single-path fit on the paired clean observation."]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--paper", type=Path, default=DEFAULT_PAPER)
    ap.add_argument("--data-dir", type=Path, default=None,
                    help="default: <paper>/gps_l1_prn09_short_delay_pilot_20261013/calibration_v2_after_acquisition_amendment")
    ap.add_argument("--out", type=Path, default=REPO / "results" / "pvt_information_test")
    ap.add_argument("--cells", nargs="*", default=list(BANKS))
    ap.add_argument("--step", type=float, default=0.002)
    ap.add_argument("--cal-boot", type=int, default=24)
    ap.add_argument("--cal-inject-boot", type=int, default=5)
    ap.add_argument("--test-boot", type=int, default=8)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    data_dir = args.data_dir or (args.paper / "gps_l1_prn09_short_delay_pilot_20261013"
                                 / "calibration_v2_after_acquisition_amendment")
    out = args.out / ("self_test" if args.self_test else "")
    out.mkdir(parents=True, exist_ok=True)
    if args.self_test:
        args.cells = args.cells[:2]
        args.cal_boot, args.cal_inject_boot, args.test_boot = 6, 1, 2
        args.step = 0.004
    jobs = [{"cell": c, "seed": SEED + i, "self_test": args.self_test, "data_dir": str(data_dir),
             "paper": str(args.paper), "step": args.step, "cal_boot": args.cal_boot,
             "cal_inject_boot": args.cal_inject_boot, "test_boot": args.test_boot}
            for i, c in enumerate(args.cells)]
    if args.workers > 1 and len(jobs) > 1:
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            results = list(pool.map(process_cell, jobs))
    else:
        results = [process_cell(j) for j in jobs]
    rows = [r for res in results for r in res["rows"]]
    calib = [res["calib"] for res in results]
    write_csv(out / "stage_a_cases.csv", rows)
    summary = summarize(rows)
    payload = {"status": "complete", "self_test": args.self_test, "iq_generated": False,
               "pvt_run": False, "design": {"delays_chip": DELAYS_CHIP, "amplitudes": AMPLITUDES,
                                            "phases_deg": PHASES_DEG},
               "settings": vars(args) | {"paper": str(args.paper), "out": str(args.out),
                                          "data_dir": str(data_dir)},
               "calibration": calib, "summary": summary}
    (out / "stage_a_summary.json").write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")
    write_summary_md(out / "STAGE_A_REPORT.md", summary, calib)
    print(f"[stage A] decision: {summary['stage_a_decision']}  -> {out}")


if __name__ == "__main__":
    main()
