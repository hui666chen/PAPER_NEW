#!/usr/bin/env python3
"""Fit two overlapping GPS C/A correlation lobes to the CUDA scan output.

The scan is referenced to the receiver's tracked prompt, so the individual
component positions have an arbitrary common offset.  Their separation does
not.  A triangular C/A main-lobe model is sufficient for this first controlled
sanity check; the next estimator can replace it with a calibrated template.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares


def read_scan(path: Path) -> tuple[np.ndarray, np.ndarray]:
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError(f"empty scan file: {path}")

    x = np.asarray([float(row["delay_chip"]) for row in rows])
    if "I_avg" in rows[0]:
        real_key, imag_key = "I_avg", "Q_avg"
    else:
        real_key, imag_key = "I", "Q"
    z = np.asarray(
        [float(row[real_key]) + 1j * float(row[imag_key]) for row in rows],
        dtype=np.complex128,
    )

    # Keep the central main-lobe region and exclude far correlation sidelobes.
    keep = (x >= -0.8) & (x <= 0.8)
    return x[keep], z[keep]


def triangular(x: np.ndarray) -> np.ndarray:
    return np.maximum(1.0 - np.abs(x), 0.0)


def fit_scan(x: np.ndarray, z: np.ndarray) -> tuple[np.ndarray, complex, complex, float]:
    def solve_linear(q: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        tau0, tau1, width = q
        design = np.column_stack(
            [
                triangular((x - tau0) / width),
                triangular((x - tau1) / width),
                np.ones_like(x),
                x,
            ]
        )
        coeff = np.linalg.lstsq(design, z, rcond=None)[0]
        error = design @ coeff - z
        residual = np.concatenate((error.real, error.imag))
        return residual, coeff

    starts = (
        (-0.15, 0.05, 0.9),
        (-0.10, 0.10, 1.0),
        (-0.25, 0.15, 1.1),
        (-0.05, 0.25, 0.8),
    )
    best = None
    for start in starts:
        result = least_squares(
            lambda q: solve_linear(q)[0],
            start,
            bounds=([-1.0, -1.0, 0.5], [1.0, 1.0, 2.0]),
            max_nfev=5000,
        )
        residual, coeff = solve_linear(result.x)
        rmse = float(np.sqrt(np.mean(residual**2)))
        if best is None or rmse < best[0]:
            best = (rmse, result.x.copy(), coeff.copy())

    assert best is not None
    rmse, parameters, coeff = best
    if parameters[0] > parameters[1]:
        parameters[[0, 1]] = parameters[[1, 0]]
        coeff[[0, 1]] = coeff[[1, 0]]
    return parameters, coeff[0], coeff[1], rmse


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("scan", type=Path)
    args = parser.parse_args()

    x, z = read_scan(args.scan)
    parameters, c0, c1, rmse = fit_scan(x, z)
    tau0, tau1, width = parameters
    separation = tau1 - tau0
    print(f"file={args.scan}")
    print(f"tau0={tau0:.5f} chip")
    print(f"tau1={tau1:.5f} chip")
    print(f"delta_tau_hat={separation:.5f} chip")
    print(f"common_lobe_width={width:.5f} chip")
    print(f"complex_amplitude_ratio={abs(c1) / abs(c0):.5f}")
    print(f"fit_rmse={rmse:.3f}")


if __name__ == "__main__":
    main()
