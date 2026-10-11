#!/usr/bin/env python3
"""Fit a dual-component scan using a measured single-component template."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares


def read_complex_scan(path: Path) -> tuple[np.ndarray, np.ndarray]:
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError(f"empty scan file: {path}")
    real_key = "I_avg" if "I_avg" in rows[0] else "I"
    imag_key = "Q_avg" if "Q_avg" in rows[0] else "Q"
    x = np.asarray([float(row["delay_chip"]) for row in rows])
    z = np.asarray(
        [float(row[real_key]) + 1j * float(row[imag_key]) for row in rows],
        dtype=np.complex128,
    )
    return x, z


def interpolate_template(
    x_template: np.ndarray, template: np.ndarray, query: np.ndarray
) -> np.ndarray:
    return np.interp(query, x_template, template.real, left=0.0, right=0.0) + 1j * np.interp(
        query, x_template, template.imag, left=0.0, right=0.0
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("template", type=Path)
    parser.add_argument("dual_scan", type=Path)
    parser.add_argument("--fit-min", type=float, default=-1.2)
    parser.add_argument("--fit-max", type=float, default=1.2)
    parser.add_argument("--min-separation", type=float, default=0.02)
    parser.add_argument("--max-separation", type=float, default=1.20)
    args = parser.parse_args()

    template_x, template_z = read_complex_scan(args.template)
    dual_x, dual_z = read_complex_scan(args.dual_scan)
    keep = (dual_x >= args.fit_min) & (dual_x <= args.fit_max)
    x = dual_x[keep]
    y = dual_z[keep]

    # With the delay-tap convention used by the CUDA scanner, a physically
    # delayed replica appears to the left of the primary component.  Fit the
    # primary position and a positive separation explicitly so that a0/a1 do
    # not get silently swapped when the positions are sorted.
    def solve_linear(parameters: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        tau_primary, separation = parameters
        tau_replica = tau_primary - separation
        basis0 = interpolate_template(template_x, template_z, x - tau_primary)
        basis1 = interpolate_template(template_x, template_z, x - tau_replica)
        design = np.column_stack([basis0, basis1, np.ones_like(x), x])
        coeff = np.linalg.lstsq(design, y, rcond=None)[0]
        error = design @ coeff - y
        return np.concatenate((error.real, error.imag)), coeff

    # Variable projection makes the complex amplitudes linear.  Rank a broad
    # coarse grid first, then refine several distinct candidates; this is much
    # less vulnerable to the local minima produced by anti-phase components.
    coarse = []
    for tau_primary in np.arange(-0.8, 0.81, 0.05):
        for separation in np.arange(args.min_separation, args.max_separation + 0.001, 0.05):
            residual, _ = solve_linear(np.asarray([tau_primary, separation]))
            coarse.append((float(np.mean(residual**2)), tau_primary, separation))

    starts = []
    for _, tau_primary, separation in sorted(coarse)[:40]:
        if all(abs(tau_primary - p) > 0.04 or abs(separation - d) > 0.04 for p, d in starts):
            starts.append((tau_primary, separation))
        if len(starts) == 12:
            break

    best = None
    for start in starts:
        result = least_squares(
            lambda parameters: solve_linear(parameters)[0],
            start,
            bounds=([-1.5, args.min_separation], [1.5, args.max_separation]),
            max_nfev=5000,
        )
        residual, coeff = solve_linear(result.x)
        rmse = float(np.sqrt(np.mean(residual**2)))
        if best is None or rmse < best[0]:
            best = (rmse, result.x.copy(), coeff.copy(), result.success, result.status)

    assert best is not None
    rmse, parameters, coeff, optimizer_success, optimizer_status = best
    tau_primary, separation = parameters
    tau_replica = tau_primary - separation
    signal_rms = float(np.sqrt(np.mean(np.concatenate((y.real, y.imag)) ** 2)))

    print(f"template={args.template}")
    print(f"dual_scan={args.dual_scan}")
    print(f"fit_range=[{args.fit_min:.2f}, {args.fit_max:.2f}] chip")
    print(f"tau_primary={tau_primary:.5f} chip")
    print(f"tau_replica={tau_replica:.5f} chip")
    print(f"delta_tau_hat={separation:.5f} chip")
    print(f"complex_amplitude_ratio={abs(coeff[1]) / abs(coeff[0]):.5f}")
    print(f"relative_phase_hat={np.angle(coeff[1] / coeff[0]):.5f} rad")
    print(f"fit_rmse={rmse:.3f}")
    print(f"fit_nrmse={rmse / signal_rms:.6f}")
    print(f"optimizer_success={str(optimizer_success).lower()}")
    print(f"optimizer_status={optimizer_status}")


if __name__ == "__main__":
    main()
