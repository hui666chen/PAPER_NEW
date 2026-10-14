#!/usr/bin/env python3
"""Deterministic known-truth audit of the historical four-column PVT fit.

This isolates the exact numerical pattern used by the PRN9/PRN30 PVT
pipelines: an unscaled inverse in the coarse profile and an unscaled
pseudoinverse in the continuous-tau objective. It does not re-score real
PVT windows, alter gates, or claim a Monte Carlo performance result.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np


PAPER = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PAPER))

from receiver_chain_response_control import ResponseSpec, load_response_model  # noqa: E402


BANKS = {
    9: PAPER / "gps_l1_prn09_high_elevation_validation/response_banks/bank_68.npz",
    30: PAPER / "gps_l1_prn30_formal_pvt/response_banks/bank_68.npz",
}
RHO = 0.34 * np.exp(1j * np.deg2rad(22.5))
TAU_OFFSETS = np.round(np.arange(-0.04, 0.04001, 0.01), 2)


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


def profile_cost(y: np.ndarray, design: np.ndarray, method: str) -> float:
    """Return summed profiled residual, with historical and stable solvers."""
    if method == "legacy_inverse":
        gram = design.conj().T @ design
        inverse = np.linalg.inv(gram)
        corr = y @ design.conj()
        explained = np.einsum(
            "ei,ij,ej->e", corr.conj(), inverse, corr, optimize=True
        ).real
        # This is the scorer's cost up to an observation-only constant.
        return float(-np.sum(explained))
    if method == "legacy_pinv":
        gram = design.conj().T @ design
        inverse = np.linalg.pinv(gram, rcond=1e-12)
        corr = y @ design.conj()
        explained = np.einsum(
            "ei,ij,ej->e", corr.conj(), inverse, corr, optimize=True
        ).real
        return float(-np.sum(explained))
    # Column scaling leaves the fitted subspace unchanged and prevents the
    # response columns (~1e7 norm) from numerically erasing the affine columns.
    scale = np.maximum(np.linalg.norm(design, axis=0), 1e-30)
    q, _ = np.linalg.qr(design / scale, mode="reduced")
    residual = y - (y @ q.conj()) @ q.T
    return float(np.sum(np.abs(residual) ** 2))


def run_case(prn: int, model, spec, separation: float, affine_fraction: float) -> dict:
    taps = np.asarray(spec.tap_grid, float)
    one = np.ones_like(taps)
    direct, reflected = model.exact(np.asarray((0.0, separation), float))
    epoch = np.arange(8)
    gains = np.exp(1j * np.linspace(-0.6, 0.7, len(epoch))) * (
        1.0 + 0.15 * np.sin(epoch)
    )
    constant = (
        affine_fraction * np.linalg.norm(direct) / np.sqrt(len(taps))
        * gains * np.exp(1j * 0.3)
    )
    linear = (
        affine_fraction * np.linalg.norm(direct) / np.linalg.norm(taps)
        * gains * np.exp(-1j * 0.2)
    )
    if affine_fraction == 0.0:
        constant = np.zeros_like(gains)
        linear = np.zeros_like(gains)
    y = (
        gains[:, None] * (direct + RHO * reflected)[None, :]
        + constant[:, None] * one[None, :]
        + linear[:, None] * taps[None, :]
    )

    tau_grid = TAU_OFFSETS
    separation_grid = np.round(
        np.arange(separation - 0.04, separation + 0.04001, 0.01), 2
    )
    candidates = {
        "legacy_coarse_inverse_4col": ("legacy_inverse", True),
        "legacy_refine_pinv_4col": ("legacy_pinv", True),
        "stable_qr_affine_4col": ("stable_qr", True),
        "stable_qr_no_affine_2col": ("stable_qr", False),
    }
    scores: dict[str, list[tuple[float, float, float]]] = {
        key: [] for key in candidates
    }
    for candidate_separation in separation_grid:
        for tau in tau_grid:
            path = model.exact(
                np.asarray((-tau, candidate_separation - tau), float)
            ).T
            design4 = np.column_stack((path, one, taps))
            design2 = path
            for key, (method, include_affine) in candidates.items():
                design = design4 if include_affine else design2
                cost = profile_cost(y, design, method)
                scores[key].append((cost, float(tau), float(candidate_separation)))

    truth_design = np.column_stack((direct, reflected, one, taps))
    gram = truth_design.conj().T @ truth_design
    eigenvalues = np.linalg.eigvalsh(gram)
    scaled_condition = float(np.linalg.cond(
        truth_design / np.maximum(np.linalg.norm(truth_design, axis=0), 1e-30)
    ))
    selected = {}
    for key, values in scores.items():
        best = min(values, key=lambda item: item[0])
        truth_cost = min(
            (cost for cost, tau, sep in values
             if abs(tau) < 1e-12 and abs(sep - separation) < 1e-12),
            default=None,
        )
        selected[key] = {
            "tau_chip": best[1],
            "separation_chip": best[2],
            "selected_offset_from_truth_chip": [best[1], best[2] - separation],
            "best_cost": best[0],
            "truth_grid_cost": truth_cost,
            "selected_cost_minus_truth_cost": (
                float(best[0] - truth_cost) if truth_cost is not None else None
            ),
        }
    return {
        "prn": prn,
        "truth": {"tau_chip": 0.0, "separation_chip": separation,
                  "rho_magnitude": float(abs(RHO)), "rho_phase_deg": 22.5},
        "affine_fraction_of_path_rms_each": affine_fraction,
        "epochs": len(epoch),
        "search": {
            "scope": "local grid only; tau +/-0.04 chip and separation +/-0.04 chip",
            "step_chip": 0.01,
            "candidates_per_solver": int(len(tau_grid) * len(separation_grid)),
        },
        "raw_gram_eigenvalues": [float(v) for v in eigenvalues],
        "raw_gram_smallest_to_largest": float(eigenvalues[0] / eigenvalues[-1]),
        "legacy_pinv_rcond": 1e-12,
        "effective_raw_gram_rank_at_truth": int(
            np.sum(eigenvalues > eigenvalues[-1] * 1e-12)
        ),
        "column_scaled_design_condition": scaled_condition,
        "solvers": selected,
    }


def main() -> None:
    assert_row_qr_projection()
    cases = []
    for prn, bank in BANKS.items():
        spec = ResponseSpec(prn=prn)
        model = load_response_model(bank, expected_spec=spec)
        for separation in (0.8, 0.1):
            for affine_fraction in (0.0, 0.05):
                cases.append(run_case(prn, model, spec, separation, affine_fraction))
    output = {
        "schema": 1,
        "status": "four_column_known_truth_noiseless_diagnostic",
        "truth_accessed": True,
        "monte_carlo_trials": 0,
        "pvt_recomputed": False,
        "cases": cases,
        "interpretation_limit": (
            "A local deterministic synthetic audit only; it does not establish "
            "the effect on historical gates or PVT. The historical pipeline is "
            "hybrid: raw np.linalg.inv in coarse candidate scoring and raw "
            "pinv(rcond=1e-12) in continuous tau refinement."
        ),
    }
    out = PAPER / "PAPER_NEW/results/shared_ratio_numerical_audit_20261014"
    out.mkdir(parents=True, exist_ok=True)
    path = out / "historical_four_column_noiseless_corrected.json"
    path.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(path)
    for case in cases:
        fields = []
        for name, result in case["solvers"].items():
            fields.append(
                f"{name}={result['tau_chip']:.2f}/{result['separation_chip']:.2f}"
            )
        print(
            f"PRN{case['prn']} sep={case['truth']['separation_chip']:.2f} "
            f"affine={case['affine_fraction_of_path_rms_each']:.2f}: "
            + ", ".join(fields)
        )


if __name__ == "__main__":
    main()
