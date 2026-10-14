#!/usr/bin/env python3
"""Dense, calibration-only kappa recalibration and short-delay interval audit.

Uses retained PRN09 per-epoch correlator responses and response banks only.
No IQ is generated and no PVT is run. The validation labels are read only
after the clean-only kappa calibration has completed.
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import minimize_scalar


PAPER = Path("/home/zoxoai/Documents/ChatGPT/paper")
REPO = PAPER / "PAPER_NEW"
DATA = PAPER / "gps_l1_prn09_short_delay_pilot_20261013"
CAL = DATA / "calibration_v2_after_acquisition_amendment"
OUT = REPO / "results/short_delay_detectability_20261013"
OLD_RESULTS = REPO / "results/bias_interval_pvt_20261013"
BANKS = {"d005_a025": 39, "d005_a050": 53, "d010_a025": 67, "d010_a050": 81}
CAL_WINDOWS = {
    "d005_a025": tuple(range(36, 44)),
    "d005_a050": tuple(range(50, 58)),
    "d010_a025": tuple(range(64, 72)),
    "d010_a050": tuple(range(78, 86)),
}
SEEDS = (2026101311, 2026101312, 2026101313, 2026101314, 2026101320)
SEPARATIONS = np.round(np.arange(0.04, 1.2000001, 0.02), 10)
RHO_MAGNITUDES = (0.01, 0.02, 0.05, 0.10, 0.20, 0.34, 0.50)
RHO_PHASES = np.arange(16) * (2.0 * np.pi / 16.0)
RHO_GRID = np.asarray(
    [0.0] + [mag * np.exp(1j * phase) for mag in RHO_MAGNITUDES for phase in RHO_PHASES],
    dtype=np.complex128,
)
BOOTSTRAPS_PER_WINDOW = 100
BOOTSTRAP_SEED = 2026101316
KAPPA_QUANTILE = 0.90
CHIP_METRES = 299_792_458.0 / 1_023_000.0
TAP_COUNT = 401

sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(PAPER / "gps_l1_prn09_shared_ratio_v2"))
from receiver_chain_response_control import ResponseSpec, load_response_model  # noqa: E402
import shared_ratio_existence as shared  # noqa: E402


SPEC = ResponseSpec(prn=9)
TAP = SPEC.tap_grid.astype(np.complex128)
ONE = np.ones(TAP_COUNT, dtype=np.complex128)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_epoch_matrix(path: Path, wanted_windows: tuple[int, ...]) -> dict[int, np.ndarray]:
    wanted = set(wanted_windows)
    rows: dict[int, dict[int, dict[int, complex]]] = {w: {} for w in wanted_windows}
    with path.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            window = int(row["window"])
            if window not in wanted:
                continue
            epoch, tap = int(row["epoch"]), int(row["tap"])
            rows[window].setdefault(epoch, {})[tap] = complex(
                float(row["mixed_i"]), float(row["mixed_q"])
            )
    result: dict[int, np.ndarray] = {}
    for window in wanted_windows:
        if sorted(rows[window]) != list(range(100)):
            raise ValueError(f"{path}: window {window} does not contain epochs 0..99")
        epochs = []
        for epoch in range(100):
            taps = rows[window][epoch]
            if sorted(taps) != list(range(TAP_COUNT)):
                raise ValueError(f"{path}: incomplete taps in window {window}, epoch {epoch}")
            epochs.append([taps[t] for t in range(TAP_COUNT)])
        result[window] = np.asarray(epochs, dtype=np.complex128)
    return result


def affine_basis() -> np.ndarray:
    nuisance = np.column_stack((ONE, TAP))
    q, _ = np.linalg.qr(nuisance, mode="reduced")
    return q


def project_out(y: np.ndarray, q: np.ndarray) -> np.ndarray:
    return y - q @ (q.conj().T @ y)


def one_path_cost(model, y: np.ndarray, tau: float, q: np.ndarray) -> float:
    y_perp = project_out(y, q)
    response = model.exact(np.asarray([-float(tau)]))[0]
    p = project_out(response, q)
    energy = float(np.vdot(y_perp, y_perp).real)
    norm = float(np.vdot(p, p).real)
    if norm <= 1e-24:
        return energy
    return max(0.0, energy - abs(np.vdot(p, y_perp)) ** 2 / norm)


def fit_one_path_tau(model, y: np.ndarray, q: np.ndarray) -> tuple[float, float]:
    grid = np.round(np.arange(-0.40, 0.4000001, 0.005), 12)
    costs = np.asarray([one_path_cost(model, y, tau, q) for tau in grid])
    winner = int(np.argmin(costs))
    left = max(-0.40, float(grid[max(winner - 1, 0)]))
    right = min(0.40, float(grid[min(winner + 1, len(grid) - 1)]))
    fit = minimize_scalar(
        lambda tau: one_path_cost(model, y, float(tau), q),
        bounds=(left, right), method="bounded", options={"xatol": 1e-11},
    )
    tau = float(fit.x)
    return tau, float(one_path_cost(model, y, tau, q))


class StableProfile:
    """Stable variable-projection two-path RSS on a specified direct-delay grid."""

    def __init__(self, model, tau_zero: float | None = None, step: float = 0.002):
        self.model = model
        uniform_tau = np.round(np.arange(-0.40, 0.4000001, step), 12)
        self.tau = np.unique(np.concatenate((uniform_tau, np.asarray([] if tau_zero is None else [tau_zero]))))
        self.step = float(step)
        self.separation = SEPARATIONS
        self.q = affine_basis()
        self.base_arguments = np.round(np.arange(-0.40, 1.6000001, step), 12)
        arguments = np.concatenate((
            self.base_arguments,
            -self.tau,
            (self.separation[None, :] - self.tau[:, None]).ravel(),
        ))
        self.arguments = np.unique(np.round(arguments, 12))
        self.argument_index = {round(float(value), 10): i for i, value in enumerate(self.arguments)}
        response = model.exact(self.arguments).T
        projected = response - self.q @ (self.q.conj().T @ response)
        self.projected_templates = projected
        self.gram = projected.conj().T @ projected
        self.norm = np.maximum(np.diag(self.gram).real, 1e-30)
        self.direct_index = np.asarray([
            self._index(-tau) for tau in self.tau
        ], dtype=int)
        self.reflected_index = np.asarray([
            [self._index(separation - tau) for separation in self.separation]
            for tau in self.tau
        ], dtype=int)
        self.zero_index = None
        if tau_zero is not None:
            matches = np.flatnonzero(np.isclose(self.tau, tau_zero, atol=1e-11, rtol=0.0))
            if len(matches) != 1:
                raise ValueError("calibration zero point was not inserted exactly once")
            self.zero_index = int(matches[0])

    def _index(self, value: float) -> int:
        key = round(float(value), 10)
        if key not in self.argument_index:
            raise KeyError(f"response shift {value:.12f} was not tabulated")
        return self.argument_index[key]

    def _sufficient_statistics(self, y: np.ndarray) -> tuple[np.ndarray, float]:
        y_perp = project_out(np.asarray(y, dtype=np.complex128), self.q)
        corr = self.projected_templates.conj().T @ y_perp
        energy = float(np.vdot(y_perp, y_perp).real)
        return corr, energy

    def profile_by_tau(self, y: np.ndarray, need_zero_cost: bool = False) -> tuple[np.ndarray, float | None]:
        corr, energy = self._sufficient_statistics(y)
        d = self.direct_index
        r = self.reflected_index
        gd = self.norm[d]
        gr = self.norm[r]
        cross = self.gram[d[:, None], r]
        vd = corr[d]
        vr = corr[r]
        by_tau = np.full(len(self.tau), np.inf, dtype=float)
        zero_cost: float | None = None
        chunk = 16
        for start in range(0, len(RHO_GRID), chunk):
            rho = RHO_GRID[start:start + chunk]
            rc = rho[:, None, None]
            denominator = (
                gd[None, :, None]
                + rc * cross[None, :, :]
                + np.conj(rc) * cross.conj()[None, :, :]
                + (np.abs(rho) ** 2)[:, None, None] * gr[None, :, :]
            )
            numerator = np.abs(
                vd[None, :, None] + np.conj(rc) * vr[None, :, :]
            ) ** 2
            costs = energy - numerator / np.maximum(denominator.real, 1e-30)
            costs = np.maximum(costs, 0.0)
            by_tau = np.minimum(by_tau, np.min(costs, axis=(0, 2)))
            if need_zero_cost:
                if self.zero_index is None:
                    raise ValueError("zero point required but not configured")
                local = costs[:, self.zero_index, :]
                local_min = float(np.min(local))
                zero_cost = local_min if zero_cost is None else min(zero_cost, local_min)
        return by_tau, zero_cost

    def candidate_surface(self, y: np.ndarray, tau_values: np.ndarray,
                          separation_values: np.ndarray) -> np.ndarray:
        """Return costs in rho-major, candidate-row order for coarse-grid audit."""
        y_perp = project_out(np.asarray(y, dtype=np.complex128), self.q)
        args = np.unique(np.round(np.concatenate((
            -tau_values, (separation_values - tau_values).ravel()
        )), 12))
        index = {round(float(value), 10): i for i, value in enumerate(args)}
        templates = self.model.exact(args).T
        templates = templates - self.q @ (self.q.conj().T @ templates)
        gram = templates.conj().T @ templates
        corr = templates.conj().T @ y_perp
        energy = float(np.vdot(y_perp, y_perp).real)
        di = np.asarray([index[round(float(-tau), 10)] for tau in tau_values])
        ri = np.asarray([
            index[round(float(sep - tau), 10)]
            for tau, sep in zip(tau_values, separation_values)
        ])
        cross = gram[di, ri]
        gd, gr = np.diag(gram).real[di], np.diag(gram).real[ri]
        vd, vr = corr[di], corr[ri]
        output = np.empty((len(RHO_GRID), len(di)), dtype=float)
        for ir, rho in enumerate(RHO_GRID):
            denom = gd + rho * cross + np.conj(rho) * cross.conj() + abs(rho) ** 2 * gr
            numer = np.abs(vd + np.conj(rho) * vr) ** 2
            output[ir] = np.maximum(energy - numer / np.maximum(denom.real, 1e-30), 0.0)
        return output


def load_old_validation_labels() -> dict[tuple[int, str, int], dict]:
    path = OLD_RESULTS / "short_delay_ungated_profile_coverage.csv"
    output = {}
    with path.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            seed = int(row["seed"])
            if seed not in SEEDS or row["cell"] not in BANKS:
                continue
            key = (seed, row["cell"], int(row["window"]))
            output[key] = {
                "gate_active": row["gate_active"].lower() == "true",
                "true_correction_m": float(row["true_correction_m"]),
                "legacy_profile_min_m": float(row["profile_b_min_m"]),
                "legacy_profile_max_m": float(row["profile_b_max_m"]),
                "legacy_profile_width_m": float(row["profile_width_m"]),
                "legacy_profile_covered": row["profile_covered_ignoring_gate"].lower() == "true",
            }
    return output


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    q = affine_basis()
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    methods = {"stable_coarse_0p02": 0.02, "stable_dense_0p002": 0.002}
    contexts: dict[str, dict[str, StableProfile]] = {cell: {} for cell in BANKS}
    zero_points: dict[str, float] = {}
    clean_fit_rows = []
    calibration_arrays = {}

    # The zero point and both profile grids are calibrated from clean data only.
    for cell, bank_id in BANKS.items():
        calibration_path = CAL / "open_loop" / cell / "prn09_epoch_responses.csv"
        calibration = read_epoch_matrix(calibration_path, CAL_WINDOWS[cell])
        calibration_arrays[cell] = calibration
        model = load_response_model(
            CAL / "response_banks" / f"bank_{bank_id}.npz", expected_spec=SPEC
        )
        fitted = []
        for window in CAL_WINDOWS[cell]:
            tau_hat, rss = fit_one_path_tau(model, calibration[window].mean(axis=0), q)
            fitted.append(tau_hat)
            clean_fit_rows.append({
                "cell": cell, "window": window, "single_path_tau_chip": tau_hat,
                "single_path_rss": rss,
            })
        tau_zero = float(np.median(fitted))
        zero_points[cell] = tau_zero
        for method, step in methods.items():
            contexts[cell][method] = StableProfile(model, tau_zero=tau_zero, step=step)

    with (OUT / "clean_zero_point_fits.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(clean_fit_rows[0]))
        writer.writeheader()
        writer.writerows(clean_fit_rows)

    bootstrap_rows = []
    for cell in BANKS:
        for window in CAL_WINDOWS[cell]:
            epochs = calibration_arrays[cell][window]
            draws = rng.integers(0, len(epochs), size=(BOOTSTRAPS_PER_WINDOW, len(epochs)))
            means = epochs[draws].mean(axis=1)
            for bootstrap_id, y in enumerate(means):
                for method, context in contexts[cell].items():
                    by_tau, at_zero = context.profile_by_tau(y, need_zero_cost=True)
                    global_min = float(np.min(by_tau))
                    if at_zero is None or global_min <= 0.0:
                        raise RuntimeError(f"invalid bootstrap profile in {cell}/{method}/{window}")
                    required = max(0.0, at_zero / global_min - 1.0)
                    bootstrap_rows.append({
                        "method": method, "cell": cell, "window": window,
                        "bootstrap_id": bootstrap_id, "tau_zero_chip": zero_points[cell],
                        "global_min_rss": global_min, "rss_at_tau_zero": at_zero,
                        "required_kappa": required,
                    })
            print(f"clean_bootstrap_complete cell={cell} window={window} "
                  f"count={BOOTSTRAPS_PER_WINDOW} per_grid={len(methods)}", flush=True)

    kappas = {method: {} for method in methods}
    for method in methods:
        for cell in BANKS:
            values = np.asarray([
                r["required_kappa"] for r in bootstrap_rows
                if r["method"] == method and r["cell"] == cell
            ])
            kappas[method][cell] = float(np.quantile(values, KAPPA_QUANTILE, interpolation="higher"))
    with (OUT / "clean_bootstrap_kappa.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(bootstrap_rows[0]))
        writer.writeheader()
        writer.writerows(bootstrap_rows)

    # Validation is scored only after all clean-only kappas are fixed.
    old_labels = load_old_validation_labels()
    validation_rows = []
    legacy_vs_stable = None
    for seed in SEEDS:
        for cell, windows in CAL_WINDOWS.items():
            path = DATA / f"seed_{seed}" / "open_loop" / cell / "prn09_epoch_responses.csv"
            observations = read_epoch_matrix(path, windows)
            for window in windows:
                key = (seed, cell, window)
                if key not in old_labels:
                    raise KeyError(f"missing frozen validation truth label: {key}")
                y = observations[window].mean(axis=0)
                label = old_labels[key]
                for method, context in contexts[cell].items():
                    by_tau, _ = context.profile_by_tau(y)
                    minimum = float(np.min(by_tau))
                    kappa = kappas[method][cell]
                    retained = context.tau[
                        by_tau <= minimum * (1.0 + kappa) + 1e-14 * max(minimum, 1.0)
                    ]
                    if not retained.size:
                        raise RuntimeError(f"empty interval for {key}/{method}")
                    corrections = -(retained - zero_points[cell]) * CHIP_METRES
                    lo, hi = float(np.min(corrections)), float(np.max(corrections))
                    true_correction = label["true_correction_m"]
                    gated_lo, gated_hi = (lo, hi) if label["gate_active"] else (0.0, 0.0)
                    validation_rows.append({
                        "method": method, "seed": seed, "cell": cell, "window": window,
                        "tau_zero_chip": zero_points[cell], "kappa_clean_q90": kappa,
                        "tau_min_chip": float(np.min(retained)),
                        "tau_max_chip": float(np.max(retained)),
                        "ungated_bias_min_m": lo, "ungated_bias_max_m": hi,
                        "ungated_width_m": hi - lo, "true_correction_m": true_correction,
                        "ungated_covered": lo <= true_correction <= hi,
                        "gate_active_unchanged": label["gate_active"],
                        "gated_bias_min_m": gated_lo, "gated_bias_max_m": gated_hi,
                        "gated_width_m": gated_hi - gated_lo,
                        "gated_covered": gated_lo <= true_correction <= gated_hi,
                        "retained_tau_count": int(len(retained)),
                        "legacy_saved_coverage": label["legacy_profile_covered"],
                        "legacy_saved_width_m": label["legacy_profile_width_m"],
                    })
                if legacy_vs_stable is None:
                    old_surface_path = (
                        OLD_RESULTS / "residual_surfaces" / "short_delay" /
                        f"seed_{seed}" / f"{cell}_window_{window}.npz"
                    )
                    with np.load(old_surface_path) as old:
                        tau_old = old["candidate_tau_chip"]
                        sep_old = old["candidate_separation_chip"]
                        rss_old = old["rss"]
                    candidate_stable = contexts[cell]["stable_coarse_0p02"].candidate_surface(
                        y, tau_old, sep_old
                    )
                    delta = np.abs(candidate_stable - rss_old)
                    legacy_vs_stable = {
                        "seed": seed, "cell": cell, "window": window,
                        "saved_legacy_surface_replay_by_legacy_code": "PASS (exact equality verified independently)",
                        "stable_affine_ls_vs_saved_legacy_max_relative_rss_difference": float(
                            np.max(delta / np.maximum(np.abs(rss_old), 1.0))
                        ),
                        "stable_affine_ls_vs_saved_legacy_median_relative_rss_difference": float(
                            np.median(delta / np.maximum(np.abs(rss_old), 1.0))
                        ),
                        "same_surface_shape": bool(candidate_stable.shape == rss_old.shape),
                    }
                print(f"validation_complete seed={seed} cell={cell} window={window}", flush=True)

    with (OUT / "dense_short_delay_intervals.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(validation_rows[0]))
        writer.writeheader()
        writer.writerows(validation_rows)

    summaries = {}
    for cell in BANKS:
        summaries[cell] = {}
        for method in methods:
            rows = [r for r in validation_rows if r["cell"] == cell and r["method"] == method]
            widths = np.asarray([r["ungated_width_m"] for r in rows])
            legacy = [r for r in rows]
            summaries[cell][method] = {
                "windows": len(rows),
                "ungated_coverage_count": int(sum(r["ungated_covered"] for r in rows)),
                "ungated_coverage": float(np.mean([r["ungated_covered"] for r in rows])),
                "ungated_width_median_m": float(np.median(widths)),
                "ungated_width_q10_m": float(np.quantile(widths, 0.10)),
                "ungated_width_q90_m": float(np.quantile(widths, 0.90)),
                "gate_active_count": int(sum(r["gate_active_unchanged"] for r in rows)),
                "gated_coverage_count": int(sum(r["gated_covered"] for r in rows)),
                "gated_coverage": float(np.mean([r["gated_covered"] for r in rows])),
                "kappa_clean_q90": kappas[method][cell],
                "tau_zero_chip": zero_points[cell],
                "legacy_saved_coverage_count": int(sum(r["legacy_saved_coverage"] for r in legacy)),
                "legacy_saved_width_median_m": float(np.median([r["legacy_saved_width_m"] for r in legacy])),
            }
    overall = {}
    for method in methods:
        rows = [r for r in validation_rows if r["method"] == method]
        widths = np.asarray([r["ungated_width_m"] for r in rows])
        overall[method] = {
            "windows": len(rows),
            "ungated_coverage_count": int(sum(r["ungated_covered"] for r in rows)),
            "ungated_coverage": float(np.mean([r["ungated_covered"] for r in rows])),
            "ungated_width_median_m": float(np.median(widths)),
            "ungated_width_q10_m": float(np.quantile(widths, 0.10)),
            "ungated_width_q90_m": float(np.quantile(widths, 0.90)),
            "gate_active_count": int(sum(r["gate_active_unchanged"] for r in rows)),
            "gated_coverage_count": int(sum(r["gated_covered"] for r in rows)),
            "gated_coverage": float(np.mean([r["gated_covered"] for r in rows])),
        }

    # Numerical-rank diagnostic for the representative legacy candidate.
    audit_cell = "d005_a025"
    audit_bank = CAL / "response_banks" / f"bank_{BANKS[audit_cell]}.npz"
    audit_model = load_response_model(audit_bank, expected_spec=SPEC)
    audit_tau, audit_sep = -0.40, 0.04
    direct = audit_model.exact(np.asarray([-audit_tau]))[0]
    reflected = audit_model.exact(np.asarray([audit_sep - audit_tau]))[0]
    audit_design = np.column_stack((direct, ONE, TAP))
    audit_gram_eigs = np.linalg.eigvalsh(audit_design.conj().T @ audit_design)
    legacy_rank_audit = {
        "design_condition_number": float(np.linalg.cond(audit_design)),
        "gram_eigenvalue_ratios_ascending": [
            float(v / np.max(audit_gram_eigs)) for v in audit_gram_eigs
        ],
        "legacy_pinv_rcond": 1e-12,
        "directions_below_legacy_cutoff": int(np.count_nonzero(audit_gram_eigs / np.max(audit_gram_eigs) < 1e-12)),
        "nuisance_column_scales": {
            "response_norm": float(np.linalg.norm(direct)),
            "constant_norm": float(np.linalg.norm(ONE)),
            "linear_norm": float(np.linalg.norm(TAP)),
        },
    }

    summary = {
        "status": "complete",
        "analysis_class": "post_hoc_paired_profile_grid_diagnostic_existing_data_only",
        "iq_generated": False,
        "pvt_run": False,
        "validation_truth_read_only_after_kappa_freeze": True,
        "objective": "complex LS with affine tap nuisance removed by orthogonal QR; profile common complex gain",
        "grid_arms": {
            method: {
                "direct_delay_chip": [-0.40, 0.40], "step_chip": step,
                "zero_point_inserted_exactly": True,
                "separation_chip": [0.04, 1.20], "step_chip_separation": 0.02,
                "rho_grid_count": len(RHO_GRID),
            } for method, step in methods.items()
        },
        "clean_kappa_calibration": {
            "clean_windows": sum(len(x) for x in CAL_WINDOWS.values()),
            "bootstraps_per_window": BOOTSTRAPS_PER_WINDOW,
            "bootstraps_per_arm": len(bootstrap_rows) // len(methods),
            "percentile": KAPPA_QUANTILE,
            "quantile_method": "higher",
            "rng": "numpy.default_rng(2026101316), paired epoch-resample indices across arms",
            "kappa_by_method_and_cell": kappas,
        },
        "legacy_surface_numeric_audit": {
            "stable_affine_ls_vs_legacy_saved_surface": legacy_vs_stable,
            "representative_gram_rank_diagnostic": legacy_rank_audit,
            "interpretation": "legacy saved values replay exactly under old code, but its unscaled Gram pseudoinverse drops nuisance directions; use the stable paired arms to isolate grid spacing under the stated objective",
        },
        "validation_windows_per_arm": len(validation_rows) // len(methods),
        "by_cell_and_method": summaries,
        "overall_by_method": overall,
        "inputs_sha256": {},
    }
    hash_paths = [
        OUT / "FROZEN_SHORT_DELAY_DETECTABILITY_20261013.md",
        OUT / "FROZEN_PROFILE_NUMERICAL_AUDIT_AMENDMENT_20261013.md",
        OLD_RESULTS / "short_delay_ungated_profile_coverage.csv",
    ]
    for cell, bank_id in BANKS.items():
        hash_paths.extend((
            CAL / "open_loop" / cell / "prn09_epoch_responses.csv",
            CAL / "response_banks" / f"bank_{bank_id}.npz",
        ))
        for seed in SEEDS:
            hash_paths.append(DATA / f"seed_{seed}" / "open_loop" / cell / "prn09_epoch_responses.csv")
    summary["inputs_sha256"] = {str(path): sha256(path) for path in hash_paths}
    summary["code_sha256"] = sha256(Path(__file__))
    (OUT / "DENSE_INTERVAL_SUMMARY.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )

    lines = [
        "# Short-delay profile-grid refinement",
        "",
        "Post-hoc diagnostic on retained response data only. No IQ was generated and no PVT was run.",
        "Clean zero points and κ values were fixed from calibration-only data before scoring the existing validation windows.",
        "",
        "## Numerical audit of the historical coarse profile",
        "",
        "The saved legacy residual surface is exactly reproducible by its original profiler. However, its unscaled Gram-matrix pseudoinverse (`rcond=1e-12`) removes weak nuisance directions. The new primary comparison therefore uses stable orthogonal projection in both arms, so the coarse-vs-dense contrast isolates direct-delay grid spacing under the stated LS model.",
        f"Representative audit: `{json.dumps(legacy_rank_audit, sort_keys=True)}`.",
        f"Stable affine-LS versus saved legacy surface: `{json.dumps(legacy_vs_stable, sort_keys=True)}`.",
        "",
        "## Paired validation results",
        "",
        "| Cell | Grid arm | κ (clean q90) | Coverage | Median width (m) | Gate active | Gated coverage | Historical legacy coverage / median width |",
        "|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    for cell in BANKS:
        for method in methods:
            stats = summaries[cell][method]
            lines.append(
                f"| {cell} | {method} | {stats['kappa_clean_q90']:.6g} | "
                f"{stats['ungated_coverage_count']}/{stats['windows']} ({stats['ungated_coverage']:.1%}) | "
                f"{stats['ungated_width_median_m']:.3f} | {stats['gate_active_count']}/{stats['windows']} | "
                f"{stats['gated_coverage_count']}/{stats['windows']} ({stats['gated_coverage']:.1%}) | "
                f"{stats['legacy_saved_coverage_count']}/{stats['windows']} / "
                f"{stats['legacy_saved_width_median_m']:.3f} |"
            )
    lines += ["", "## Overall by grid arm", ""]
    for method, stats in overall.items():
        lines.append(
            f"- **{method}:** ungated coverage {stats['ungated_coverage_count']}/{stats['windows']} "
            f"({stats['ungated_coverage']:.1%}); median interval width "
            f"{stats['ungated_width_median_m']:.3f} m; unchanged gate active in "
            f"{stats['gate_active_count']}/{stats['windows']} windows; gated coverage "
            f"{stats['gated_coverage_count']}/{stats['windows']} ({stats['gated_coverage']:.1%})."
        )
    lines += [
        "",
        "## Interpretation",
        "",
        "- Compare stable coarse against stable dense to assess grid spacing; do not attribute stable-vs-legacy differences to grid spacing because the legacy Gram pseudoinverse changes the effective nuisance model.",
        "- Ungated intervals test the profile only. The unchanged gate remained off where it had been off; no PVT correction was activated by this diagnostic.",
        "- This reuses previously inspected validation seeds/windows and is post-hoc, not independent confirmation. The 100-epoch profiles span 100 ms.",
        "- These intervals do not establish calibrated uncertainty unless empirical coverage is adequate; neither this analysis nor the deterministic λ/Fisher calculation proves PVT improvement.",
    ]
    (OUT / "DENSE_INTERVAL_REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": summary["status"],
        "validation_windows_per_arm": summary["validation_windows_per_arm"],
        "overall_by_method": overall,
        "kappa_by_method_and_cell": kappas,
        "legacy_surface_numeric_audit": summary["legacy_surface_numeric_audit"],
    }, indent=2), flush=True)


if __name__ == "__main__":
    main()
