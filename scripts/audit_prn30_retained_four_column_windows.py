#!/usr/bin/env python3
"""Re-score retained PRN30 observations with historical and stable profiles.

This is a numerical audit only. It does not change or rerun the PVT solver.
The per-window estimate and confidence gate are recomputed for three profile
solvers while keeping the response banks, observations, grid, Sobol settings,
and activation threshold fixed.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import minimize_scalar


PAPER = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER / "gps_l1_prn30_formal_pvt"))

from b_joint_blind_audit import TAU  # noqa: E402
from crossfit_basin_confidence import CrossfitBasinScorer  # noqa: E402
from gps_l1_prn30_formal_pvt.estimate_formal_windows import load_windows  # noqa: E402
from monte_carlo_resolution_baseline import build_dictionary  # noqa: E402
from receiver_chain_response_control import ResponseSpec, load_response_model  # noqa: E402


ROOT = PAPER / "gps_l1_prn30_formal_pvt"
OUT = PAPER / "PAPER_NEW/results/shared_ratio_numerical_audit_20261014"
SPEC = ResponseSpec(prn=30)
SCENES = (
    ("clean", tuple(range(60, 76)), 68),
    ("amp_0p01", tuple(range(79, 95)), 87),
    ("amp_0p05", tuple(range(100, 116)), 108),
    ("amp_0p10", tuple(range(121, 137)), 129),
)
GATE_THRESHOLD = 0.20
STABLE_INVERSE_QR_CHECKS = {}


def assert_row_qr_projection() -> None:
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
    error = np.linalg.norm(projected - reference) / max(np.linalg.norm(reference), 1e-30)
    if error > 1e-12:
        raise AssertionError(f"complex QR projection check failed: {error}")


def stable_profile_inverses(dictionary: np.ndarray, indices: np.ndarray) -> np.ndarray:
    """Map column-scaled stable Gram inverses back to raw dictionary coords."""
    design = np.stack([dictionary[:, row] for row in indices])
    scales = np.maximum(np.linalg.norm(design, axis=1), 1e-30)
    scaled = design / scales[:, None, :]
    gram = np.einsum("nmc,nmd->ncd", scaled.conj(), scaled, optimize=True)
    identity = np.broadcast_to(np.eye(gram.shape[-1]), gram.shape)
    inverse_scaled = np.linalg.solve(gram, identity)
    return inverse_scaled / scales[:, :, None] / scales[:, None, :]


def validate_profile_inverse(dictionary, indices, inverse, observation):
    """Compare representative scaled-Gram scores against direct complex QR."""
    selected = sorted(set((0, len(indices) // 2, len(indices) - 1)))
    errors = []
    for candidate in selected:
        design = np.stack([dictionary[:, col] for col in indices[candidate]], axis=1)
        correlation = observation @ design.conj()
        explained_gram = np.vdot(
            correlation, inverse[candidate] @ correlation
        ).real
        scale = np.maximum(np.linalg.norm(design, axis=0), 1e-30)
        q, _ = np.linalg.qr(design / scale, mode="reduced")
        explained_qr = np.linalg.norm(observation @ q.conj()) ** 2
        relative_error = abs(explained_gram - explained_qr) / max(
            abs(explained_qr), 1.0
        )
        errors.append(float(relative_error))
    if max(errors) > 1e-9:
        raise AssertionError(f"scaled-Gram/QR score mismatch: {errors}")
    return {"candidate_indices": selected, "relative_errors": errors,
            "maximum_relative_error": max(errors)}


def refine_objective(model, y, weights, separation, tau, mode):
    response = model.exact(np.asarray((-tau, separation - tau), float)).T
    if mode == "stable_affine":
        design = np.column_stack((response, np.ones(401), SPEC.tap_grid))
    elif mode == "stable_no_affine":
        design = response
    else:
        design = np.column_stack((response, np.ones(401), SPEC.tap_grid))
        gram = design.conj().T @ design
        inverse = np.linalg.pinv(gram, rcond=1e-12)
        correlations = y @ design.conj()
        explained = np.einsum(
            "ei,ij,ej->e", correlations.conj(), inverse, correlations,
            optimize=True,
        ).real
        return -float(np.sum(weights * explained)), float(np.linalg.cond(gram))

    scale = np.maximum(np.linalg.norm(design, axis=0), 1e-30)
    q, _ = np.linalg.qr(design / scale, mode="reduced")
    residual = y - (y @ q.conj()) @ q.T
    return float(np.sum(weights[:, None] * np.abs(residual) ** 2)), float(
        np.linalg.cond(design / scale)
    )


def score_scene(scene: str, windows: tuple[int, ...], bank_window: int,
                scene_index: int) -> list[dict]:
    source = ROOT / "final_sign_persistence/open_loop" / scene / "prn30_epoch_responses.csv"
    bank_path = ROOT / "response_banks" / f"bank_{bank_window}.npz"
    observations = load_windows(source, windows)
    model = load_response_model(bank_path, expected_spec=SPEC)
    _, dictionary, indices, inverse_a, _ = build_dictionary(model, SPEC)
    indices_c = indices[:, :2]
    inverse_b = stable_profile_inverses(dictionary, indices)
    inverse_c = stable_profile_inverses(dictionary, indices_c)
    validation_y = observations[windows[0]][0][0]
    STABLE_INVERSE_QR_CHECKS[scene] = {
        "with_affine": validate_profile_inverse(
            dictionary, indices, inverse_b, validation_y
        ),
        "without_affine": validate_profile_inverse(
            dictionary, indices_c, inverse_c, validation_y
        ),
    }

    solvers = {
        "A_historical_hybrid": (indices, inverse_a, "historical"),
        "B_stable_with_affine": (indices, inverse_b, "stable_affine"),
        "C_stable_without_affine": (indices_c, inverse_c, "stable_no_affine"),
    }
    records = []
    for window in windows:
        y, weights, _ = observations[window]
        for mode_name, (mode_indices, mode_inverse, refine_mode) in solvers.items():
            scorer = CrossfitBasinScorer(
                dictionary, mode_indices, mode_inverse,
                competitors=8, sobol_points=4096, sobol_scrambles=2,
            )
            corr = y @ dictionary.conj()
            score = scorer.score(
                corr, weights,
                seed=2026092900 + scene_index * 100 + window,
            )
            anchor = int(score["full_anchor"])
            coarse_tau = float(score["full_tau"])
            separation = float(scorer.candidate_separation[anchor])
            low = max(float(TAU[0]), coarse_tau - 0.01)
            high = min(float(TAU[-1]), coarse_tau + 0.01)
            coarse_value, coarse_condition = refine_objective(
                model, y, weights, separation, coarse_tau, refine_mode
            )
            fit = minimize_scalar(
                lambda tau: refine_objective(
                    model, y, weights, separation, float(tau), refine_mode
                )[0],
                bounds=(low, high), method="bounded",
                options={"xatol": 1e-6, "maxiter": 80},
            )
            refined_tau = float(fit.x)
            refined_value, refined_condition = refine_objective(
                model, y, weights, separation, refined_tau, refine_mode
            )
            objective_ok = bool(
                refined_value <= coarse_value + 1e-7 * max(1.0, abs(coarse_value))
            )
            active = bool(
                score["P_cross"] <= GATE_THRESHOLD + 1e-12
                and fit.success and objective_ok
                and not score["consistency_failure"]
            )
            records.append({
                "scene": scene, "window": int(window), "bank_window": bank_window,
                "solver": mode_name,
                "P_cross": float(score["P_cross"]),
                "consistency_failure": bool(score["consistency_failure"]),
                "cross_active": active,
                "coarse_tau_chip": coarse_tau,
                "separation_chip": separation,
                "refined_tau_chip": refined_tau,
                "correction_m": -refined_tau * (299_792_458.0 / 1.023e6),
                "coarse_score": coarse_value,
                "refined_score": refined_value,
                "objective_nonincrease": objective_ok,
                "coarse_condition": coarse_condition,
                "refined_condition": refined_condition,
                "full_tau_boundary": bool(score["full_tau_boundary"]),
                "P_A_to_B": float(score["P_A_to_B"]),
                "P_B_to_A": float(score["P_B_to_A"]),
            })
        if (window - windows[0] + 1) % 4 == 0:
            print(f"{scene}: completed {window - windows[0] + 1}/{len(windows)} windows",
                  flush=True)
    return records


def main() -> None:
    assert_row_qr_projection()
    rows = []
    for scene_index, (scene, windows, bank_window) in enumerate(SCENES):
        rows.extend(score_scene(scene, windows, bank_window, scene_index))
    csv_path = OUT / "prn30_retained_four_column_rescore_corrected.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    by_scene_solver = {}
    for scene, _, _ in SCENES:
        for solver in ("A_historical_hybrid", "B_stable_with_affine",
                       "C_stable_without_affine"):
            subset = [r for r in rows if r["scene"] == scene and r["solver"] == solver]
            by_scene_solver[f"{scene}/{solver}"] = {
                "windows": len(subset),
                "cross_active": sum(r["cross_active"] for r in subset),
                "consistency_failures": sum(r["consistency_failure"] for r in subset),
                "tau_mean_chip": float(np.mean([r["refined_tau_chip"] for r in subset])),
                "tau_std_chip": float(np.std([r["refined_tau_chip"] for r in subset], ddof=1)),
                "separation_mean_chip": float(np.mean([r["separation_chip"] for r in subset])),
                "separation_boundary_count": sum(
                    r["separation_chip"] in (0.04, 1.2) for r in subset
                ),
            }
    summary = {
        "schema": 1,
        "status": "retained_PRN30_four_column_rescore_no_PVT",
        "pvt_recomputed": False,
        "input_scope": "retained PRN30 final-sign-persistence 100x401 epoch responses",
        "frozen": {
            "scenes_windows_banks": [
                {"scene": s, "windows": [w[0], w[-1]], "bank_window": b}
                for s, w, b in SCENES
            ],
            "gate_threshold": GATE_THRESHOLD,
            "tau_grid_from_existing_pipeline": [float(TAU[0]), float(TAU[-1]), 0.01],
            "separation_grid": [0.04, 1.2, 0.01],
            "sobol_points": 4096,
            "sobol_scrambles": 2,
        },
        "solver_definitions": {
            "A_historical_hybrid": "raw-inverse coarse scorer; raw pinv(rcond=1e-12) continuous tau refinement",
            "B_stable_with_affine": "column-scaled stable profile for coarse scoring and continuous refinement, with constant/linear nuisance columns",
            "C_stable_without_affine": "column-scaled stable profile, omitting constant/linear nuisance columns",
        },
        "stable_inverse_vs_direct_qr_checks": STABLE_INVERSE_QR_CHECKS,
        "summary": by_scene_solver,
        "caution": (
            "Exploratory numerical re-score of one retained PRN30 batch. "
            "The frozen persistence aggregation and PVT were not recomputed; "
            "the 10-run PRN9 validation's full epoch responses were deleted, "
            "so this does not adjudicate the PRN9 claim."
        ),
    }
    json_path = OUT / "prn30_retained_four_column_rescore_summary_corrected.json"
    json_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
