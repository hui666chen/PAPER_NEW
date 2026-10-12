#!/usr/bin/env python3
"""Profiled M0/M1/M2 CRB grid with explicit C/N0 scaling.

This is the public, reproducible front end for the CRB/identifiability map.
The supplied noise covariance must be the independently estimated 401-tap
covariance at ``--reference-cn0``.  No simulator truth is read.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from correlator_baselines import ChainBankTemplate


DERIVATIVE_STEP = 1.0e-4


def _values(text: str, cast=float) -> list:
    return [cast(item) for item in text.split(",") if item]


def _response(model: ChainBankTemplate, delay: float) -> tuple[np.ndarray, np.ndarray]:
    values = np.stack([
        model.response(delay - DERIVATIVE_STEP),
        model.response(delay),
        model.response(delay + DERIVATIVE_STEP),
    ])
    return values[1], (values[2] - values[0]) / (2.0 * DERIVATIVE_STEP)


def _residualized_information(whiten, nuisance: np.ndarray, derivatives: np.ndarray) -> np.ndarray:
    x = whiten(nuisance)
    d = whiten(derivatives)
    coefficient = np.linalg.lstsq(x, d, rcond=1.0e-12)[0]
    residual = d - x @ coefficient
    return 2.0 * np.real(residual.conj().T @ residual)


def _finalize(fisher: np.ndarray, names: tuple[str, ...]) -> dict[str, object]:
    fisher = (fisher + fisher.T) / 2.0
    eigenvalues = np.linalg.eigvalsh(fisher)
    scale = max(float(eigenvalues[-1]), 1.0e-300)
    rank = int(np.sum(eigenvalues > scale * 1.0e-10))
    condition = (
        float(eigenvalues[-1] / eigenvalues[0])
        if eigenvalues[0] > 0.0 else float("inf")
    )
    covariance = None
    if rank == len(names):
        covariance = np.linalg.inv(fisher)
    output: dict[str, object] = {
        "rank": rank,
        "dimension": len(names),
        "condition": condition,
        "fisher_eigenvalues": [float(value) for value in eigenvalues],
    }
    for name in ("tau", "delta", "ratio", "phi0", "omega"):
        if name not in names or covariance is None:
            output[f"{name}_sd_crb"] = None
            continue
        index = names.index(name)
        output[f"{name}_sd_crb"] = float(np.sqrt(max(covariance[index, index], 0.0)))
    return output


def evaluate(
    model: ChainBankTemplate,
    whiten,
    model_name: str,
    delay: float,
    ratio: float,
    phase: float,
    epochs: int,
) -> dict[str, object]:
    r0, rp0 = _response(model, 0.0)
    r1, rp1 = _response(model, delay)
    background = np.column_stack((np.ones(model.tap_count), model.tap_grid))
    c = ratio * np.exp(1j * phase)

    if model_name == "M0":
        nuisance = np.column_stack((r0, r1, background))
        derivatives = np.column_stack((-rp0 - c * rp1, c * rp1))
        fisher = epochs * _residualized_information(whiten, nuisance, derivatives)
        names = ("tau", "delta")
    else:
        times = (np.arange(epochs) - (epochs - 1) / 2.0) * 0.001
        fisher = None
        for time_s in times:
            angle = phase
            ck = ratio * np.exp(1j * angle)
            h = r0 + ck * r1
            nuisance = np.column_stack((h, background))
            columns = [
                -rp0 - ck * rp1,
                ck * rp1,
                np.exp(1j * angle) * r1,
                1j * ck * r1,
            ]
            if model_name == "M2":
                columns.append(1j * time_s * ck * r1)
            part = _residualized_information(
                whiten, nuisance, np.column_stack(columns)
            )
            fisher = part if fisher is None else fisher + part
        names = (
            ("tau", "delta", "ratio", "phi0")
            if model_name == "M1"
            else ("tau", "delta", "ratio", "phi0", "omega")
        )

    output = _finalize(fisher, names)
    output.update({
        "model": model_name,
        "epochs": epochs,
        "delay_chip": delay,
        "amplitude_ratio": ratio,
        "phi0_rad": phase,
    })
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bank", type=Path, required=True)
    parser.add_argument("--noise-covariance", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--reference-cn0", type=float, required=True)
    parser.add_argument("--delays", default="0.05,0.10,0.20,0.40,0.60,0.80,1.00")
    parser.add_argument("--ratios", default="0.10,0.20,0.34,0.50,0.70")
    parser.add_argument("--phases-deg", default="0,90,180,270")
    parser.add_argument("--cn0", default="35,40,45,50")
    parser.add_argument("--epochs", default="20,50,100")
    parser.add_argument("--models", default="M0,M1,M2")
    args = parser.parse_args()

    model = ChainBankTemplate(args.bank)
    with np.load(args.noise_covariance, allow_pickle=False) as data:
        eigenvalues = np.asarray(data["eigenvalues"], dtype=float)
        eigenvectors = np.asarray(data["eigenvectors"], dtype=np.complex128)
    if eigenvalues.shape != (model.tap_count,) or eigenvectors.shape != (
        model.tap_count, model.tap_count
    ):
        raise ValueError("noise covariance dimensions do not match the response bank")
    if np.any(eigenvalues <= 0.0):
        raise ValueError("noise covariance eigenvalues must be positive")

    delays = _values(args.delays)
    ratios = _values(args.ratios)
    phases = [np.deg2rad(value) for value in _values(args.phases_deg)]
    cn0_values = _values(args.cn0)
    epochs_values = _values(args.epochs, int)
    models = [value for value in args.models.split(",") if value]
    invalid = sorted(set(models) - {"M0", "M1", "M2"})
    if invalid:
        raise ValueError(f"unknown models: {invalid}")

    results = []
    for cn0 in cn0_values:
        variance_scale = 10.0 ** ((args.reference_cn0 - cn0) / 10.0)

        def whiten(matrix: np.ndarray) -> np.ndarray:
            return (eigenvectors.conj().T @ matrix) / np.sqrt(
                eigenvalues * variance_scale
            )[:, None]

        for epochs in epochs_values:
            for ratio in ratios:
                for delay in delays:
                    for phase in phases:
                        for name in models:
                            row = evaluate(
                                model, whiten, name, delay, ratio, phase, epochs
                            )
                            row.update({
                                "cn0_db_hz": cn0,
                                "reference_cn0_db_hz": args.reference_cn0,
                                "noise_variance_scale": variance_scale,
                            })
                            results.append(row)
        print(f"cn0={cn0:g} complete", flush=True)

    payload = {
        "schema": 1,
        "truth_accessed": False,
        "bank": str(args.bank),
        "noise_covariance": str(args.noise_covariance),
        "derivative_step_chip": DERIVATIVE_STEP,
        "grid": {
            "delays_chip": delays,
            "amplitude_ratios": ratios,
            "phase_deg": _values(args.phases_deg),
            "cn0_db_hz": cn0_values,
            "epochs": epochs_values,
            "models": models,
        },
        "identifiability_diagnostic": {
            "raw_outputs": ["tau_sd_crb", "delta_sd_crb", "rank", "condition"],
            "optional_boundary_line": "delta_sd_crb <= 0.02 chip",
            "warning": "the line is a diagnostic visualization threshold, not a universal physical limit",
        },
        "results": results,
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"rows": len(results), "output": str(args.output_json)}))


if __name__ == "__main__":
    main()
