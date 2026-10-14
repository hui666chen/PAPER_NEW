#!/usr/bin/env python3
"""Known-truth numerical audit of the shared-ratio profile solver.

This is a diagnostic only. It does not modify the frozen estimator, gate,
correction schedules, or PVT outputs. A is the historical Gram-pseudoinverse
implementation; B profiles the same affine nuisance subspace using a QR
projection; C profiles only the composite path response with no affine terms.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np

PAPER = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER / "gps_l1_prn09_shared_ratio_v2"))

import shared_ratio_existence as old  # noqa: E402
from b_joint_blind_audit import SEPARATION, TAU  # noqa: E402
from gps_l1_prn30_formal_pvt.estimate_formal_windows import load_windows  # noqa: E402
from monte_carlo_resolution_baseline import build_dictionary  # noqa: E402
from receiver_chain_response_control import (  # noqa: E402
    ResponseSpec,
    load_response_model,
)


CALIBRATION = {
    9: {
        "bank": PAPER / "gps_l1_prn09_comparable_v3/response_banks/bank_68.npz",
        "responses": PAPER / "gps_l1_prn09_comparable_v3/calibration/gpu_reference/window_68/prn09_epoch_responses.csv",
    },
    30: {
        "bank": PAPER / "gps_l1_prn30_shared_ratio_v2/response_banks_25m/bank_68.npz",
        "responses": PAPER / "gps_l1_prn30_shared_ratio_v2/calibration_25m/gpu_reference/window_68/prn30_epoch_responses.csv",
    },
}


def assert_row_qr_projection() -> None:
    """Guard the complex row-vector projection against conjugation mistakes."""
    rng = np.random.default_rng(20261014)
    design = rng.normal(size=(31, 4)) + 1j * rng.normal(size=(31, 4))
    design[:, 2:] *= 1e-4
    observations = rng.normal(size=(5, 31)) + 1j * rng.normal(size=(5, 31))
    scaled = design / np.linalg.norm(design, axis=0)
    q, _ = np.linalg.qr(scaled, mode="reduced")
    projected = (observations @ q.conj()) @ q.T
    reference = (scaled @ np.linalg.lstsq(
        scaled, observations.T, rcond=None
    )[0]).T
    relative_error = np.linalg.norm(projected - reference) / max(
        np.linalg.norm(reference), 1e-30
    )
    if relative_error > 1e-12:
        raise AssertionError(f"complex QR projection check failed: {relative_error}")


def empirical_affine_ratios(prn: int, model, tap_grid: np.ndarray) -> dict:
    """Estimate clean-window affine scale relative to the fitted direct gain.

    These are empirical projection scales, not claims that the physical
    receiver has a literal DC/linear bias. The fit is column-scaled before LS.
    """
    data = load_windows(CALIBRATION[prn]["responses"], (68,))[68]
    values, _, _ = data
    direct = model.exact([0.0])[0]
    design = np.column_stack((direct, np.ones(len(tap_grid)), tap_grid))
    scales = np.linalg.norm(design, axis=0)
    scaled = design / scales
    ratios = []
    for y in values:
        coefficients = np.linalg.lstsq(scaled, y, rcond=None)[0] / scales
        gain, constant, linear = coefficients
        if abs(gain) > 1e-20:
            ratios.append((constant / gain, linear / gain))
    if not ratios:
        raise RuntimeError(f"PRN{prn}: could not estimate empirical affine scale")
    ratios = np.asarray(ratios, complex)
    # Componentwise complex median is deterministic and robust to a few bad epochs.
    c_ratio = complex(np.median(ratios[:, 0].real), np.median(ratios[:, 0].imag))
    d_ratio = complex(np.median(ratios[:, 1].real), np.median(ratios[:, 1].imag))
    return {
        "c_over_gain": c_ratio,
        "d_over_gain": d_ratio,
        "coefficient_ratios": ratios,
        "n_epochs": int(len(ratios)),
        "median_abs_c_over_gain": float(np.median(np.abs(ratios[:, 0]))),
        "median_abs_d_over_gain": float(np.median(np.abs(ratios[:, 1]))),
        "location_constant_rms_fraction": float(
            abs(c_ratio) * np.sqrt(len(tap_grid)) / np.linalg.norm(direct)
        ),
        "location_linear_rms_fraction": float(
            abs(d_ratio) * np.linalg.norm(tap_grid) / np.linalg.norm(direct)
        ),
        "median_abs_constant_rms_fraction": float(
            np.median(np.abs(ratios[:, 0])) * np.sqrt(len(tap_grid))
            / np.linalg.norm(direct)
        ),
        "median_abs_linear_rms_fraction": float(
            np.median(np.abs(ratios[:, 1])) * np.linalg.norm(tap_grid)
            / np.linalg.norm(direct)
        ),
        "constant_phase_concentration": float(
            abs(np.mean(ratios[:, 0])) / max(np.mean(np.abs(ratios[:, 0])), 1e-30)
        ),
        "linear_phase_concentration": float(
            abs(np.mean(ratios[:, 1])) / max(np.mean(np.abs(ratios[:, 1])), 1e-30)
        ),
    }


def stable_profile(
    values: np.ndarray,
    weights: np.ndarray,
    dictionary: np.ndarray,
    candidate_rows: np.ndarray,
    candidate_tau: np.ndarray,
    candidate_separation: np.ndarray,
    rho_values: np.ndarray,
    include_affine: bool,
) -> dict:
    """Variable projection with direct residual evaluation and no Gram pinv."""
    y = np.asarray(values, complex)
    w = np.asarray(weights, float)
    if include_affine:
        nuisance = dictionary[:, -2:]
        q, _ = np.linalg.qr(nuisance, mode="reduced")
        yp = y - (y @ q.conj()) @ q.T
    else:
        q = None
        yp = y

    direct = dictionary[:, candidate_rows[:, 0]].T
    reflected = dictionary[:, candidate_rows[:, 1]].T
    best: dict | None = None
    chunk_size = 128
    for rho in rho_values:
        for first in range(0, len(candidate_rows), chunk_size):
            last = min(first + chunk_size, len(candidate_rows))
            h = direct[first:last] + rho * reflected[first:last]
            if q is not None:
                hp = h - (h @ q.conj()) @ q.T
            else:
                hp = h
            denominator = np.einsum("cm,cm->c", hp.conj(), hp).real
            numerator = hp.conj() @ yp.T
            coefficients = numerator / np.maximum(denominator[:, None], 1e-300)
            residual = yp[None, :, :] - coefficients[:, :, None] * hp[:, None, :]
            score = np.einsum("e,cem->c", w, np.abs(residual) ** 2, optimize=True)
            local = int(np.argmin(score))
            if best is None or score[local] < best["rss"]:
                index = first + local
                best = {
                    "tau_chip": float(candidate_tau[index]),
                    "separation_chip": float(candidate_separation[index]),
                    "rho_real": float(rho.real),
                    "rho_imag": float(rho.imag),
                    "rho_magnitude": float(abs(rho)),
                    "rho_phase_deg": float(np.degrees(np.angle(rho))) if abs(rho) else 0.0,
                    "rss": float(score[local]),
                }
    if best is None:
        raise RuntimeError("empty stable profile")
    return best


def stable_profile_fast(
    values: np.ndarray,
    weights: np.ndarray,
    dictionary: np.ndarray,
    candidate_rows: np.ndarray,
    candidate_tau: np.ndarray,
    candidate_separation: np.ndarray,
    rho_values: np.ndarray,
    include_affine: bool,
) -> dict:
    """QR-projected variable projection using stable sufficient statistics.

    Used only for noisy trials whose residual energy is well above roundoff.
    The noiseless tests continue to use direct residual vectors.
    """
    y = np.asarray(values, complex)
    w = np.asarray(weights, float)
    if include_affine:
        nuisance = dictionary[:, -2:]
        q, _ = np.linalg.qr(nuisance, mode="reduced")
        y = y - (y @ q.conj()) @ q.T
        p_all = dictionary[:, candidate_rows[:, 0]].T
        r_all = dictionary[:, candidate_rows[:, 1]].T
        p_all = p_all - (p_all @ q.conj()) @ q.T
        r_all = r_all - (r_all @ q.conj()) @ q.T
    else:
        p_all = dictionary[:, candidate_rows[:, 0]].T
        r_all = dictionary[:, candidate_rows[:, 1]].T

    p_corr = (p_all.conj() @ y.T).T
    r_corr = (r_all.conj() @ y.T).T
    p_norm = np.einsum("cm,cm->c", p_all.conj(), p_all).real
    r_norm = np.einsum("cm,cm->c", r_all.conj(), r_all).real
    cross = np.einsum("cm,cm->c", p_all.conj(), r_all)
    energy = np.einsum("e,em,em->", w, y.conj(), y).real
    rhos = np.asarray(rho_values, complex)
    denominator = (
        p_norm[:, None]
        + 2.0 * np.real(cross[:, None] * rhos[None, :])
        + r_norm[:, None] * np.abs(rhos[None, :]) ** 2
    )
    numerator = (
        p_corr[:, :, None]
        + r_corr[:, :, None] * np.conj(rhos)[None, None, :]
    )
    explained = np.sum(w[:, None, None] * np.abs(numerator) ** 2, axis=0)
    score = energy - explained / np.maximum(denominator, 1e-300)
    candidate_index, rho_index = np.unravel_index(np.argmin(score), score.shape)
    rho = complex(rhos[rho_index])
    return {
        "tau_chip": float(candidate_tau[candidate_index]),
        "separation_chip": float(candidate_separation[candidate_index]),
        "rho_real": float(rho.real),
        "rho_imag": float(rho.imag),
        "rho_magnitude": float(abs(rho)),
        "rho_phase_deg": float(np.degrees(np.angle(rho))) if abs(rho) else 0.0,
        "rss": float(score[candidate_index, rho_index]),
    }


def empirical_clean_residual_pool(prn: int, exclude_window: int = 68) -> dict:
    """Cross-window empirical clean residuals, decomposed by stable LS.

    Residual vectors and nuisance coefficients are kept paired within their
    original 5-epoch blocks for the later block-resampling diagnostic.
    They are empirical receiver/model perturbations, not asserted to be pure
    thermal noise.
    """
    spec = ResponseSpec(prn=prn)
    base_bank = CALIBRATION[prn]["bank"]
    base_csv = CALIBRATION[prn]["responses"]
    bank_dir = base_bank.parent
    response_parent = base_csv.parents[1]
    pools = []
    diagnostics = []
    for window in (68, 108, 129):
        if window == exclude_window:
            continue
        model = load_response_model(
            bank_dir / f"bank_{window}.npz", expected_spec=spec
        )
        response_path = response_parent / f"window_{window}" / base_csv.name
        values, _, _ = load_windows(response_path, (window,))[window]
        direct = model.exact([0.0])[0]
        design = np.column_stack((direct, np.ones(401), spec.tap_grid.astype(float)))
        scale = np.linalg.norm(design, axis=0)
        scaled = design / scale
        residuals = []
        affine = []
        gain_values = []
        relative_norms = []
        for y in values:
            beta = np.linalg.lstsq(scaled, y, rcond=None)[0] / scale
            gain, constant, linear = beta
            if abs(gain) <= 1e-20:
                continue
            residual = y - design @ beta
            residuals.append(residual / gain)
            affine.append((constant / gain, linear / gain))
            gain_values.append(gain)
            relative_norms.append(np.linalg.norm(residual) / max(np.linalg.norm(design @ beta), 1e-30))
        residuals = np.asarray(residuals, complex)
        affine = np.asarray(affine, complex)
        if residuals.shape != (100, 401):
            raise RuntimeError(f"PRN{prn} window {window}: incomplete residual pool")
        pools.append({"window": window, "residual": residuals, "affine": affine})
        diagnostics.append({
            "window": window,
            "epochs": int(len(residuals)),
            "median_relative_clean_residual": float(np.median(relative_norms)),
            "p90_relative_clean_residual": float(np.quantile(relative_norms, 0.9)),
            "median_abs_constant_over_gain": float(np.median(np.abs(affine[:, 0]))),
            "median_abs_linear_over_gain": float(np.median(np.abs(affine[:, 1]))),
        })
    return {"pools": pools, "diagnostics": diagnostics}


def draw_empirical_blocks(pool: dict, rng: np.random.Generator,
                          epoch_count: int = 100, block_length: int = 5):
    blocks = []
    for item in pool["pools"]:
        for first in range(0, 100 - block_length + 1, block_length):
            last = first + block_length
            blocks.append((item["residual"][first:last], item["affine"][first:last]))
    needed = int(np.ceil(epoch_count / block_length))
    selected = rng.integers(0, len(blocks), size=needed)
    residual = np.concatenate([blocks[i][0] for i in selected], axis=0)[:epoch_count]
    affine = np.concatenate([blocks[i][1] for i in selected], axis=0)[:epoch_count]
    return residual, affine


def run_noisy_trials(prn: int, separation_true: float, rho_true: complex,
                     trials: int, seed: int, pool: dict) -> dict:
    """Semi-parametric 100-epoch trials from held-out clean residual blocks."""
    spec = ResponseSpec(prn=prn)
    model = load_response_model(CALIBRATION[prn]["bank"], expected_spec=spec)
    delay_grid, dictionary, _, _, _ = build_dictionary(model, spec)
    observations = load_windows(CALIBRATION[prn]["responses"], (68,))[68]
    _, weights, magnitude = observations
    weights = np.asarray(weights, float)
    gains = np.exp(1j * np.zeros(100))
    tau_true = 0.0
    direct_index = int(round((-tau_true - delay_grid[0]) / 0.01))
    reflected_index = int(round((separation_true - tau_true - delay_grid[0]) / 0.01))
    direct = dictionary[:, direct_index]
    reflected = dictionary[:, reflected_index]
    rows, candidate_tau, candidate_separation = make_candidates(
        dictionary, delay_grid, tau_true, separation_true
    )
    rho_values = old.make_rho_grid()
    records = []
    for trial in range(int(trials)):
        rng = np.random.default_rng(seed + trial)
        residual, affine = draw_empirical_blocks(pool, rng)
        gains = np.exp(1j * rng.uniform(-np.pi, np.pi, 100))
        y = np.asarray([
            gains[k] * (direct + rho_true * reflected + affine[k, 0]
                        + affine[k, 1] * spec.tap_grid + residual[k])
            for k in range(100)
        ])
        a = old.profile_shared_ratio(
            y, weights, dictionary, delay_grid, rows, candidate_tau,
            candidate_separation, rho_values,
        )
        b = stable_profile_fast(
            y, weights, dictionary, rows, candidate_tau,
            candidate_separation, rho_values, True,
        )
        c = stable_profile_fast(
            y, weights, dictionary, rows, candidate_tau,
            candidate_separation, rho_values, False,
        )
        row = {"trial": trial}
        for label, fit in (("A", a), ("B", b), ("C", c)):
            fitted_rho = complex(fit["rho_real"], fit["rho_imag"])
            row[label] = {
                "tau_chip": float(fit["tau_chip"]),
                "separation_chip": float(fit["separation_chip"]),
                "rho_magnitude": float(fit["rho_magnitude"]),
                "rho_phase_deg": float(fit["rho_phase_deg"]),
                "rss": float(fit["rss"]),
                "tau_error_chip": float(fit["tau_chip"] - tau_true),
                "separation_error_chip": float(fit["separation_chip"] - separation_true),
                "rho_abs_error": float(abs(fitted_rho - rho_true)),
            }
        records.append(row)
        if (trial + 1) % 25 == 0:
            print(f"PRN{prn} delta={separation_true:g} empirical-residual MC {trial+1}/{trials}", flush=True)

    summary = {}
    for label in ("A", "B", "C"):
        for metric in ("tau_error_chip", "separation_error_chip", "rho_abs_error"):
            errors = np.asarray([row[label][metric] for row in records], float)
            summary.setdefault(label, {})[metric] = {
                "bias": float(errors.mean()),
                "rmse": float(np.sqrt(np.mean(errors ** 2))),
                "median_abs": float(np.median(np.abs(errors))),
            }
    return {
        "prn": prn,
        "separation_true_chip": separation_true,
        "rho_true": {"real": float(rho_true.real), "imag": float(rho_true.imag),
                     "magnitude": float(abs(rho_true)),
                     "phase_deg": float(np.degrees(np.angle(rho_true)))},
        "trials": int(trials),
        "epochs_per_trial": 100,
        "noise_model": "5-epoch block resampling of paired held-out clean-window stable-LS nuisance coefficients and residual vectors, normalized by direct gain",
        "noise_pool_diagnostics": pool["diagnostics"],
        "summary": summary,
        "records": records,
    }


def make_candidates(dictionary: np.ndarray, delay_grid: np.ndarray,
                    tau_true: float, separation_true: float):
    rows, tau, separation = old.build_candidates(dictionary, delay_grid)
    keep = ((np.abs(tau - tau_true) <= 0.10 + 1e-12)
            & (np.abs(separation - separation_true) <= 0.12 + 1e-12))
    return rows[keep], tau[keep], separation[keep]


def run_case(prn: int, separation_true: float, affine_mode: str,
             rho_true: complex, seed: int) -> dict:
    spec = ResponseSpec(prn=prn)
    model = load_response_model(CALIBRATION[prn]["bank"], expected_spec=spec)
    delay_grid, dictionary, _, _, _ = build_dictionary(model, spec)
    affine_scale = empirical_affine_ratios(prn, model, spec.tap_grid.astype(float))

    tau_true = 0.0
    direct_index = int(round((-tau_true - delay_grid[0]) / 0.01))
    reflected_index = int(round((separation_true - tau_true - delay_grid[0]) / 0.01))
    direct = dictionary[:, direct_index]
    reflected = dictionary[:, reflected_index]

    rng = np.random.default_rng(seed)
    epoch_count = 8
    gains = (0.9 + 0.2 * rng.random(epoch_count)) * np.exp(
        1j * rng.uniform(-np.pi, np.pi, epoch_count)
    )
    if affine_mode == "none":
        c_by_epoch = np.zeros(epoch_count, complex)
        d_by_epoch = np.zeros(epoch_count, complex)
    elif affine_mode == "coherent_median":
        c_by_epoch = np.full(epoch_count, affine_scale["c_over_gain"], complex)
        d_by_epoch = np.full(epoch_count, affine_scale["d_over_gain"], complex)
    elif affine_mode == "empirical_resampled":
        ratios = affine_scale["coefficient_ratios"]
        picked = rng.choice(len(ratios), size=epoch_count, replace=True)
        c_by_epoch = ratios[picked, 0]
        d_by_epoch = ratios[picked, 1]
    else:
        raise ValueError(f"unknown affine mode: {affine_mode}")
    y = np.asarray([
        gain * (direct + rho_true * reflected + c_by_epoch[k]
                + d_by_epoch[k] * spec.tap_grid)
        for k, gain in enumerate(gains)
    ])
    weights = np.ones(epoch_count, float)
    rows, tau, separation = make_candidates(
        dictionary, delay_grid, tau_true, separation_true
    )
    rho_values = old.make_rho_grid()

    # A intentionally retains the historical implementation and arithmetic.
    result_a = old.profile_shared_ratio(
        y, weights, dictionary, delay_grid, rows, tau, separation, rho_values
    )
    result_b = stable_profile(
        y, weights, dictionary, rows, tau, separation, rho_values, True
    )
    result_c = stable_profile(
        y, weights, dictionary, rows, tau, separation, rho_values, False
    )

    def error(item):
        fitted_rho = complex(item["rho_real"], item["rho_imag"])
        return {
            "tau_error_chip": float(item["tau_chip"] - tau_true),
            "separation_error_chip": float(item["separation_chip"] - separation_true),
            "rho_abs_error": float(abs(fitted_rho - rho_true)),
        }

    # The old Gram matrix's nominal eigenvalue ratio is recorded only as a
    # diagnostic; roundoff can make its smallest eigenvalue slightly negative.
    representative_h = direct + rho_true * reflected
    old_design = np.column_stack((representative_h, np.ones(401), spec.tap_grid))
    gram = old_design.conj().T @ old_design
    eigenvalues = np.linalg.eigvalsh((gram + gram.conj().T) / 2)
    effective_pinv_rank = int(np.count_nonzero(
        eigenvalues > max(float(eigenvalues[-1]), 0.0) * 1e-12
    ))
    return {
        "prn": prn,
        "tau_true_chip": tau_true,
        "separation_true_chip": separation_true,
        "rho_true": {"real": float(rho_true.real), "imag": float(rho_true.imag),
                     "magnitude": float(abs(rho_true)),
                     "phase_deg": float(np.degrees(np.angle(rho_true)))},
        "affine_injection_mode": affine_mode,
        "empirical_clean_projection_scale": {
            "c_over_gain_abs_median": affine_scale["median_abs_c_over_gain"],
            "d_over_gain_abs_median": affine_scale["median_abs_d_over_gain"],
            "coherent_location_constant_rms_fraction_of_direct": affine_scale["location_constant_rms_fraction"],
            "coherent_location_linear_rms_fraction_of_direct": affine_scale["location_linear_rms_fraction"],
            "median_absolute_constant_rms_fraction_of_direct": affine_scale["median_abs_constant_rms_fraction"],
            "median_absolute_linear_rms_fraction_of_direct": affine_scale["median_abs_linear_rms_fraction"],
            "constant_phase_concentration": affine_scale["constant_phase_concentration"],
            "linear_phase_concentration": affine_scale["linear_phase_concentration"],
            "epochs": affine_scale["n_epochs"],
        },
        "candidate_count": int(len(rows)),
        "gram_eigenvalues": [float(v) for v in eigenvalues],
        "gram_pinv_rcond": 1e-12,
        "effective_pinv_rank_at_truth": effective_pinv_rank,
        "A_old_pinv": {**result_a, "error": error(result_a)},
        "B_stable_with_affine": {**result_b, "error": error(result_b)},
        "C_stable_without_affine": {**result_c, "error": error(result_c)},
    }


def main() -> None:
    assert_row_qr_projection()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prns", type=int, nargs="+", default=[9, 30])
    parser.add_argument("--separations", type=float, nargs="+", default=[0.8, 0.1])
    parser.add_argument("--seed", type=int, default=20261014)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--noise-trials", type=int, default=0,
                        help="empirical clean-residual block-bootstrap trials per PRN/delay")
    args = parser.parse_args()
    rho = 0.34 * np.exp(1j * np.deg2rad(22.5))
    rows = []
    for prn in args.prns:
        for separation in args.separations:
            for case_index, affine_mode in enumerate(
                ("none", "coherent_median", "empirical_resampled")
            ):
                rows.append(run_case(
                    prn, separation, affine_mode, rho,
                    args.seed + prn * 100 + int(separation * 1000) + 10 * case_index,
                ))
    payload = {
        "schema": 1,
        "status": "synthetic_known_truth_noiseless_only",
        "truth_accessed": True,
        "monte_carlo_noise_trials": 0,
        "search_scope": "local diagnostic neighborhood around injected tau and separation; not the full PVT estimator search",
        "rho_grid_truth_on_grid": True,
        "cases": rows,
    }
    if args.noise_trials:
        pools = {
            prn: empirical_clean_residual_pool(prn, exclude_window=68)
            for prn in args.prns
        }
        monte_carlo = []
        for prn in args.prns:
            for separation in args.separations:
                monte_carlo.append(run_noisy_trials(
                    prn, separation, rho,
                    args.noise_trials,
                    args.seed + prn * 1000 + int(separation * 10000),
                    pools[prn],
                ))
        payload["status"] = "synthetic_known_truth_with_empirical_clean_residual_trials"
        payload["noise_trials_per_prn_delay"] = int(args.noise_trials)
        payload["monte_carlo"] = monte_carlo
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": payload["status"],
        "output": str(args.output.resolve()) if args.output else None,
        "search_scope": payload["search_scope"],
        "cases": [
            {
                "prn": row["prn"],
                "separation_true_chip": row["separation_true_chip"],
                "affine_injection_mode": row["affine_injection_mode"],
                "affine_rms_fractions": [
                    row["empirical_clean_projection_scale"]["coherent_location_constant_rms_fraction_of_direct"],
                    row["empirical_clean_projection_scale"]["coherent_location_linear_rms_fraction_of_direct"],
                ],
                "median_absolute_affine_rms_fractions": [
                    row["empirical_clean_projection_scale"]["median_absolute_constant_rms_fraction_of_direct"],
                    row["empirical_clean_projection_scale"]["median_absolute_linear_rms_fraction_of_direct"],
                ],
                "effective_pinv_rank_at_truth": row["effective_pinv_rank_at_truth"],
                "fit": {
                    name: [item["tau_chip"], item["separation_chip"],
                           item["rho_magnitude"], item["rho_phase_deg"]]
                    for name, item in (
                        ("A", row["A_old_pinv"]),
                        ("B", row["B_stable_with_affine"]),
                        ("C", row["C_stable_without_affine"]),
                    )
                },
                "errors": {
                    name: item["error"]
                    for name, item in (
                        ("A", row["A_old_pinv"]),
                        ("B", row["B_stable_with_affine"]),
                        ("C", row["C_stable_without_affine"]),
                    )
                },
                "rss": {
                    name: item["rss"]
                    for name, item in (
                        ("A", row["A_old_pinv"]),
                        ("B", row["B_stable_with_affine"]),
                        ("C", row["C_stable_without_affine"]),
                    )
                },
            }
            for row in rows
        ],
        "monte_carlo_summary": [
            {
                "prn": row["prn"],
                "separation_true_chip": row["separation_true_chip"],
                "trials": row["trials"],
                "summary": row["summary"],
                "noise_pool_diagnostics": row["noise_pool_diagnostics"],
            }
            for row in payload.get("monte_carlo", [])
        ],
    }, indent=2))


if __name__ == "__main__":
    main()
