#!/usr/bin/env python3
"""Fit a two-component delay model to every continuous CUDA scan window.

The GPU block writes one 401-tap scan per row group.  This script reuses a
single-component, phase-aligned scan as the measured template and performs a
variable-projection fit for the primary delay and positive replica separation
in every window.  The nonlinear covariance of the separation is reported as
the per-measurement uncertainty used by the Kalman comparison.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares


def load_scan(path: Path) -> tuple[np.ndarray, np.ndarray]:
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        real_key = "I_avg" if "I_avg" in (reader.fieldnames or []) else "I"
        imag_key = "Q_avg" if "Q_avg" in (reader.fieldnames or []) else "Q"
        rows = [
            row for row in reader
            if all(str(row.get(key, "") or "").strip() for key in ("delay_chip", real_key, imag_key))
        ]
    if not rows:
        raise ValueError(f"empty scan file: {path}")
    x = np.asarray([float(row["delay_chip"]) for row in rows], dtype=np.float64)
    z = np.asarray(
        [float(row[real_key]) + 1j * float(row[imag_key]) for row in rows],
        dtype=np.complex128,
    )
    order = np.argsort(x)
    return x[order], z[order]


def load_windows(path: Path) -> list[tuple[int, float, np.ndarray, np.ndarray]]:
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        rows = [
            row for row in reader
            if all(str(row.get(key, "") or "").strip() for key in
                   ("window", "time_ms", "delay_chip", "I_avg", "Q_avg", "mag_avg"))
        ]
    required = {"window", "time_ms", "delay_chip", "I_avg", "Q_avg", "mag_avg"}
    if not rows:
        raise ValueError(f"empty window scan file: {path}")
    missing = required.difference(rows[0])
    if missing:
        raise ValueError(f"{path} is missing columns: {sorted(missing)}")

    grouped: dict[int, list[dict[str, str]]] = {}
    for row in rows:
        grouped.setdefault(int(row["window"]), []).append(row)
    result = []
    for window, group in sorted(grouped.items()):
        group.sort(key=lambda row: float(row["delay_chip"]))
        time_values = np.asarray([float(row["time_ms"]) for row in group])
        x = np.asarray([float(row["delay_chip"]) for row in group])
        z = np.asarray(
            [float(row["I_avg"]) + 1j * float(row["Q_avg"]) for row in group],
            dtype=np.complex128,
        )
        if len(x) < 8 or not np.all(np.isfinite(x)) or not np.all(np.isfinite(z)):
            continue
        result.append((window, float(np.mean(time_values)), x, z))
    if not result:
        raise ValueError(f"no usable windows in {path}")
    return result


def interpolate_complex(x_template: np.ndarray, template: np.ndarray, query: np.ndarray) -> np.ndarray:
    return np.interp(query, x_template, template.real, left=0.0, right=0.0) + 1j * np.interp(
        query, x_template, template.imag, left=0.0, right=0.0
    )


def fit_window(
    x_template: np.ndarray,
    template: np.ndarray,
    x_scan: np.ndarray,
    y_scan: np.ndarray,
    previous: np.ndarray | None,
    fit_min: float,
    fit_max: float,
    min_separation: float,
    max_separation: float,
) -> dict[str, float | bool]:
    keep = (x_scan >= fit_min) & (x_scan <= fit_max)
    x = x_scan[keep]
    y = y_scan[keep]
    if len(x) < 8:
        raise ValueError("too few samples in fit interval")

    def solve_linear(parameters: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        tau_primary, separation = parameters
        tau_replica = tau_primary - separation
        basis_primary = interpolate_complex(x_template, template, x - tau_primary)
        basis_replica = interpolate_complex(x_template, template, x - tau_replica)
        design = np.column_stack([basis_primary, basis_replica, np.ones_like(x), x])
        coefficients = np.linalg.lstsq(design, y, rcond=None)[0]
        residual_complex = design @ coefficients - y
        return np.concatenate((residual_complex.real, residual_complex.imag)), coefficients

    starts: list[tuple[float, float]] = []
    if previous is not None and np.all(np.isfinite(previous)):
        starts.append((float(previous[0]), float(np.clip(previous[1], min_separation, max_separation))))
    starts.extend(
        [
            (-0.20, 0.10),
            (0.00, 0.20),
            (0.20, 0.30),
            (0.00, 0.45),
            (-0.35, 0.50),
        ]
    )

    best: tuple[float, np.ndarray, np.ndarray, object] | None = None
    for start in starts:
        clipped_start = np.asarray(
            [np.clip(start[0], -1.5, 1.5), np.clip(start[1], min_separation, max_separation)],
            dtype=np.float64,
        )
        result = least_squares(
            lambda parameters: solve_linear(parameters)[0],
            clipped_start,
            bounds=([-1.5, min_separation], [1.5, max_separation]),
            max_nfev=800,
        )
        residual, coefficients = solve_linear(result.x)
        rmse = float(np.sqrt(np.mean(residual**2)))
        if best is None or rmse < best[0]:
            best = (rmse, result.x.copy(), coefficients.copy(), result)
    assert best is not None
    rmse, parameters, coefficients, optimizer = best
    tau_primary, separation = parameters
    signal_rms = float(np.sqrt(np.mean(np.concatenate((y.real, y.imag)) ** 2)))

    # The least-squares Jacobian is with respect to the two nonlinear
    # parameters after the complex amplitudes and linear nuisance terms have
    # been projected out.  This gives a local covariance for separation.
    dof = max(len(optimizer.fun) - 6, 1)
    residual_variance = float(np.sum(optimizer.fun**2) / dof)
    normal = optimizer.jac.T @ optimizer.jac
    covariance = np.linalg.pinv(normal, rcond=1.0e-12) * residual_variance
    sigma = float(math.sqrt(max(float(covariance[1, 1]), 0.0)))
    # A zero residual can occur for a synthetic/template-identical window; a
    # small floor keeps the likelihood numerically meaningful.
    sigma = max(sigma, 1.0e-4)
    ratio = abs(coefficients[1]) / max(abs(coefficients[0]), 1.0e-12)
    phase = np.angle(coefficients[1] / coefficients[0]) if abs(coefficients[0]) > 1.0e-12 else float("nan")
    return {
        "tau_primary": float(tau_primary),
        "tau_replica": float(tau_primary - separation),
        "delta_tau_hat": float(separation),
        "sigma_chip": sigma,
        "measurement_variance_chip2": sigma * sigma,
        "amplitude_ratio_hat": float(ratio),
        "relative_phase_hat_rad": float(phase),
        "fit_rmse": rmse,
        "fit_nrmse": rmse / max(signal_rms, 1.0e-12),
        "fit_success": bool(optimizer.success),
    }


def fit_prn(
    prn: int,
    kind: str,
    root: Path,
    template_root: Path,
    args: argparse.Namespace,
) -> list[dict[str, object]]:
    scan_path = root / f"{kind}_prn{prn}" / f"gps_prn{prn}_cuda_delay_scan_windows.csv"
    template_path = template_root / f"template_prn{prn}" / f"gps_prn{prn}_cuda_delay_scan_phase_aligned_avg.csv"
    if not scan_path.exists():
        raise FileNotFoundError(f"missing continuous scan: {scan_path}")
    if not template_path.exists():
        raise FileNotFoundError(f"missing template scan: {template_path}")
    template_x, template_z = load_scan(template_path)
    output = []
    previous = None
    for window, time_ms, x, z in load_windows(scan_path):
        fit = fit_window(
            template_x,
            template_z,
            x,
            z,
            previous,
            args.fit_min,
            args.fit_max,
            args.min_separation,
            args.max_separation,
        )
        previous = np.asarray([fit["tau_primary"], fit["delta_tau_hat"]], dtype=np.float64)
        output.append({"scenario": kind, "prn": prn, "window": window, "time_ms": time_ms, **fit})
    return output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True, help="multisat experiment root")
    parser.add_argument("--template-root", type=Path, default=None)
    parser.add_argument("--kind", choices=("spoof", "multipath"), required=True)
    parser.add_argument("--prns", nargs="+", type=int, default=[3, 4, 7, 8, 9, 14, 16, 27])
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--fit-min", type=float, default=-1.2)
    parser.add_argument("--fit-max", type=float, default=1.2)
    parser.add_argument("--min-separation", type=float, default=0.04)
    parser.add_argument("--max-separation", type=float, default=0.80)
    args = parser.parse_args()
    template_root = args.template_root or args.root
    output_path = args.output or args.root / f"{args.kind}_window_estimates.csv"

    rows: list[dict[str, object]] = []
    for prn in args.prns:
        rows.extend(fit_prn(prn, args.kind, args.root, template_root, args))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "scenario", "prn", "window", "time_ms", "delta_tau_hat", "sigma_chip",
        "measurement_variance_chip2", "tau_primary", "tau_replica", "amplitude_ratio_hat",
        "relative_phase_hat_rad", "fit_rmse", "fit_nrmse", "fit_success",
    ]
    with output_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    print(f"scenario={args.kind}")
    print(f"prns={','.join(str(prn) for prn in args.prns)}")
    print(f"windows={len(rows) // max(len(args.prns), 1)}")
    print(f"output={output_path}")


if __name__ == "__main__":
    main()
