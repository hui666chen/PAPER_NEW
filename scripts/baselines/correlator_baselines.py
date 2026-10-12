#!/usr/bin/env python3
"""Common-observation baselines for the 401-point GPS L1 C/A scan.

This module deliberately keeps the proposed receiver method out of the
baseline implementations.  It supplies three classical comparators:

* EML/NEML: a coherent early-minus-late discriminator with 0.1-chip spacing;
* a normalized double-delta discriminator with frozen inner/outer spacings;
* a two-component MEDLL-like variable-projection fit.

The MEDLL fit can use either an ideal triangular C/A lobe or the calibrated
receiver-chain response bank.  Running both with identical grids, nuisance
terms, and optimizer settings is the controlled ``ideal template`` versus
``receiver-chain template`` ablation required by the paper protocol.

The code consumes the scanner's epoch-response CSV and never reads simulator
truth.  Truth may only be joined later by a separate scoring script.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import json
from pathlib import Path
from typing import Iterable

import numpy as np
from scipy.optimize import minimize


CHIP_METRES = 299_792_458.0 / 1.023e6


@dataclass(frozen=True)
class BaselineConfig:
    """Frozen observation-level settings.

    The delay signs follow the current receiver convention: a component with
    primary delay ``tau`` is represented by the receiver response at ``-tau``;
    a reflected component with separation ``delta`` is at ``delta - tau``.
    """

    eml_spacing_chip: float = 0.10
    dd_inner_spacing_chip: float = 0.10
    dd_outer_spacing_chip: float = 0.50
    discriminator_min_chip: float = -0.50
    discriminator_max_chip: float = 0.50
    discriminator_step_chip: float = 0.002
    fit_min_chip: float = -1.50
    fit_max_chip: float = 1.50
    medll_tau_min_chip: float = -0.50
    medll_tau_max_chip: float = 0.50
    medll_tau_step_chip: float = 0.02
    medll_separation_min_chip: float = 0.02
    medll_separation_max_chip: float = 1.00
    medll_separation_step_chip: float = 0.04
    medll_refine_starts: int = 8
    medll_max_nfev: int = 400


@dataclass(frozen=True)
class EpochScan:
    window: int
    epoch: int
    bank_window: int
    delay_chip: np.ndarray
    mixed: np.ndarray
    metadata: dict[str, str]


def _first(row: dict[str, str], names: Iterable[str]) -> str:
    for name in names:
        if name in row and row[name] != "":
            return row[name]
    raise KeyError(f"none of {tuple(names)} found in CSV row")


def read_epoch_scans(
    path: Path, expected_taps: int = 401, default_bank_window: int | None = None
) -> list[EpochScan]:
    """Read and group one scanner epoch-response CSV.

    Accepted IQ column names are the current ``mixed_i/mixed_q`` names and
    the older ``I_avg/Q_avg`` or ``I/Q`` names.  No truth or scene parameter
    is required.
    """

    groups: dict[tuple[int, int], list[dict[str, str]]] = {}
    with path.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            window = int(float(row["window"]))
            epoch = int(float(row["epoch"]))
            groups.setdefault((window, epoch), []).append(row)
    scans: list[EpochScan] = []
    for (window, epoch), rows in sorted(groups.items()):
        rows.sort(key=lambda row: float(row["delay_chip"]))
        x = np.asarray([float(row["delay_chip"]) for row in rows], dtype=float)
        real = np.asarray(
            [float(_first(row, ("mixed_i", "I_avg", "I"))) for row in rows],
            dtype=float,
        )
        imag = np.asarray(
            [float(_first(row, ("mixed_q", "Q_avg", "Q"))) for row in rows],
            dtype=float,
        )
        if expected_taps and len(x) != expected_taps:
            raise ValueError(
                f"{path}: window={window}, epoch={epoch} has {len(x)} taps, "
                f"expected {expected_taps}"
            )
        if len(x) < 5 or np.any(np.diff(x) <= 0):
            raise ValueError(f"{path}: non-increasing or too-short delay grid")
        bank_value = rows[0].get("bank_window", "")
        bank_window = int(float(bank_value)) if bank_value else (
            window if default_bank_window is None else default_bank_window
        )
        metadata = {
            key: value for key, value in rows[0].items()
            if key not in {"mixed_i", "mixed_q", "I_avg", "Q_avg", "I", "Q"}
        }
        scans.append(EpochScan(window, epoch, bank_window, x, real + 1j * imag, metadata))
    if not scans:
        raise ValueError(f"empty epoch-response CSV: {path}")
    return scans


def _interp(x: np.ndarray, y: np.ndarray, query: np.ndarray) -> np.ndarray:
    return np.interp(query, x, y, left=0.0, right=0.0)


def _nearest_root(grid: np.ndarray, values: np.ndarray) -> tuple[float, bool]:
    """Return the root nearest zero and whether the search hit a boundary."""

    finite = np.isfinite(values)
    candidates: list[float] = []
    for index in range(len(grid) - 1):
        left, right = values[index], values[index + 1]
        if not finite[index] or not finite[index + 1]:
            continue
        if left == 0.0:
            candidates.append(float(grid[index]))
        if left * right < 0.0:
            weight = abs(left) / (abs(left) + abs(right))
            candidates.append(float(grid[index] + weight * (grid[index + 1] - grid[index])))
    if finite[-1] and values[-1] == 0.0:
        candidates.append(float(grid[-1]))
    if candidates:
        estimate = min(candidates, key=abs)
    else:
        estimate = float(grid[np.nanargmin(np.abs(values))])
    boundary = bool(np.isclose(estimate, grid[0]) or np.isclose(estimate, grid[-1]))
    return estimate, boundary


def _discriminator(
    x: np.ndarray,
    z: np.ndarray,
    spacing: float,
    normalize_spacing: bool,
    cfg: BaselineConfig,
) -> tuple[float, bool, float]:
    grid = np.arange(
        cfg.discriminator_min_chip,
        cfg.discriminator_max_chip + cfg.discriminator_step_chip * 0.5,
        cfg.discriminator_step_chip,
    )
    power = np.abs(z) ** 2
    early = _interp(x, power, grid - spacing / 2.0)
    late = _interp(x, power, grid + spacing / 2.0)
    value = early - late
    if normalize_spacing:
        outer = spacing
        value = value / outer
    estimate, boundary = _nearest_root(grid, value)
    residual = float(np.interp(estimate, grid, value))
    return estimate, boundary, residual


def eml_estimate(scan: EpochScan, cfg: BaselineConfig) -> dict[str, float | bool]:
    estimate, boundary, residual = _discriminator(
        scan.delay_chip, scan.mixed, cfg.eml_spacing_chip, False, cfg
    )
    return {
        "tau_chip": estimate,
        "boundary": boundary,
        "discriminator_residual": residual,
        "spacing_chip": cfg.eml_spacing_chip,
    }


def double_delta_estimate(
    scan: EpochScan, cfg: BaselineConfig
) -> dict[str, float | bool]:
    """Normalized double-delta discriminator.

    The two early-minus-late differences are divided by their spacings before
    subtraction.  This removes the first-order spacing scale from the
    discriminator and makes the exact comparator reproducible; it is not
    presented as a new discriminator design.
    """

    x, z = scan.delay_chip, scan.mixed
    grid = np.arange(
        cfg.discriminator_min_chip,
        cfg.discriminator_max_chip + cfg.discriminator_step_chip * 0.5,
        cfg.discriminator_step_chip,
    )
    power = np.abs(z) ** 2
    di = cfg.dd_inner_spacing_chip
    do = cfg.dd_outer_spacing_chip
    inner = (_interp(x, power, grid - di / 2.0) - _interp(x, power, grid + di / 2.0)) / di
    outer = (_interp(x, power, grid - do / 2.0) - _interp(x, power, grid + do / 2.0)) / do
    value = inner - outer
    estimate, boundary = _nearest_root(grid, value)
    return {
        "tau_chip": estimate,
        "boundary": boundary,
        "discriminator_residual": float(np.interp(estimate, grid, value)),
        "inner_spacing_chip": di,
        "outer_spacing_chip": do,
    }


class IdealTemplate:
    """Ideal triangular GPS C/A correlation template."""

    def __init__(self, x: np.ndarray, center_chip: float):
        self.center_chip = float(center_chip)
        self._x = np.asarray(x, dtype=float)

    def component(self, x: np.ndarray, delay_chip: float) -> np.ndarray:
        return np.maximum(1.0 - np.abs(x - (self.center_chip + delay_chip)), 0.0).astype(
            np.complex128
        )


class ChainBankTemplate:
    """Standalone loader for the public receiver-chain response bank format."""

    def __init__(self, path: Path):
        with np.load(path, allow_pickle=False) as data:
            self.bank = np.asarray(data["bank"], dtype=np.complex128)
            self.kmin = int(np.asarray(data["kmin"]).reshape(-1)[0])
            self.kmax = int(np.asarray(data["kmax"]).reshape(-1)[0])
            raw = np.asarray(data["metadata"]).reshape(-1)[0]
            self.metadata = json.loads(str(raw))
        spec = self.metadata.get("spec", {})
        self.sample_rate_hz = float(spec["sample_rate_hz"])
        self.if_hz = float(spec["if_hz"])
        self.chip_rate_hz = float(spec["chip_rate_hz"])
        self.half = int(spec.get("delay_filter_half", 128))
        self.beta = float(spec.get("delay_filter_beta", 10.0))
        self.chip_per_sample = self.chip_rate_hz / self.sample_rate_hz
        self.omega_per_sample = 2.0 * np.pi * self.if_hz / self.sample_rate_hz
        self.tap_first_chip = float(spec.get("tap_first_chip", -2.0))
        self.tap_step_chip = float(spec.get("tap_step_chip", 0.01))
        self.tap_count = int(spec.get("tap_count", self.bank.shape[1]))
        self.tap_grid = self.tap_first_chip + self.tap_step_chip * np.arange(self.tap_count)
        self._response_cache: dict[float, np.ndarray] = {}
        if self.bank.shape != (self.kmax - self.kmin + 1, self.tap_count):
            raise ValueError(f"{path}: inconsistent bank shape {self.bank.shape}")

    def response(self, delay_chip: float) -> np.ndarray:
        cache_key = float(delay_chip)
        cached = self._response_cache.get(cache_key)
        if cached is not None:
            return cached
        delay_samples = float(delay_chip) / self.chip_per_sample
        integer = int(np.floor(delay_samples))
        fraction = delay_samples - integer
        offsets = np.arange(-self.half, self.half + 1)
        indices = integer + offsets - self.kmin
        if indices.min() < 0 or indices.max() >= len(self.bank):
            raise ValueError(f"delay {delay_chip} chip outside response-bank support")
        kernel = np.sinc(offsets - fraction) * np.kaiser(2 * self.half + 1, self.beta)
        kernel /= kernel.sum()
        value = kernel @ self.bank[indices]
        response = value * np.exp(1j * self.omega_per_sample * delay_samples)
        self._response_cache[cache_key] = response
        return response

    def component(self, x: np.ndarray, delay_chip: float) -> np.ndarray:
        # This is the convention used by the existing shared-ratio estimator:
        # primary response(-tau), reflected response(delta-tau).
        full = self.response(-float(delay_chip))
        return _interp(self.tap_grid, full.real, x) + 1j * _interp(
            self.tap_grid, full.imag, x
        )

    def center_chip(self) -> float:
        return float(self.tap_grid[np.argmax(np.abs(self.response(0.0)))])


def _linear_fit(template, x: np.ndarray, y: np.ndarray, tau: float, separation: float):
    primary = template.component(x, tau)
    reflected = template.component(x, tau - separation)
    design = np.column_stack((primary, reflected, np.ones_like(x), x))
    scales = np.linalg.norm(design, axis=0)
    scales = np.maximum(scales, 1.0e-30)
    scaled = design / scales
    coefficient_scaled, *_ = np.linalg.lstsq(scaled, y, rcond=1.0e-12)
    coefficient = coefficient_scaled / scales
    residual = y - design @ coefficient
    return float(np.vdot(residual, residual).real), coefficient, residual


def _batch_linear_costs(
    template, x: np.ndarray, y: np.ndarray,
    tau_grid: np.ndarray, separation: float,
) -> np.ndarray:
    """Profile all direct-delay candidates for one separation at once."""

    primary = np.stack(
        [template.component(x, float(tau)) for tau in tau_grid], axis=0
    )
    reflected = np.stack(
        [template.component(x, float(tau - separation)) for tau in tau_grid], axis=0
    )
    constant = np.broadcast_to(np.ones_like(x), primary.shape)
    linear = np.broadcast_to(x, primary.shape)
    design = np.stack((primary, reflected, constant, linear), axis=2)
    scales = np.sqrt(np.einsum("nmc,nmc->nc", design.conj(), design).real)
    scales = np.maximum(scales, 1.0e-30)
    scaled = design / scales[:, None, :]
    gram = np.einsum("nmc,nmd->ncd", scaled.conj(), scaled)
    rhs = np.einsum("nmc,m->nc", scaled.conj(), y)
    try:
        scaled_coefficients = np.linalg.solve(gram, rhs[..., None])[..., 0]
    except np.linalg.LinAlgError:
        scaled_coefficients = np.einsum(
            "ncd,nd->nc", np.linalg.pinv(gram, rcond=1.0e-12), rhs
        )
    eigenvalues = np.linalg.eigvalsh((gram + gram.conj().transpose(0, 2, 1)) * 0.5)
    unstable = (~np.isfinite(eigenvalues[:, 0])) | (
        eigenvalues[:, 0] < np.maximum(eigenvalues[:, -1], 1.0e-30) * 1.0e-12
    )
    coefficients = scaled_coefficients / scales
    residual = y[None, :] - np.einsum("nmc,nc->nm", design, coefficients)
    costs = np.sum(np.abs(residual) ** 2, axis=1).real
    for index in np.flatnonzero(unstable):
        costs[index] = _linear_fit(
            template, x, y, float(tau_grid[index]), float(separation)
        )[0]
    return costs


def medll_estimate(
    scan: EpochScan, template, cfg: BaselineConfig
) -> dict[str, float | bool]:
    """Fit a two-component MEDLL-like model with variable projection."""

    keep = (scan.delay_chip >= cfg.fit_min_chip) & (scan.delay_chip <= cfg.fit_max_chip)
    x, y = scan.delay_chip[keep], scan.mixed[keep]
    tau_grid = np.arange(
        cfg.medll_tau_min_chip,
        cfg.medll_tau_max_chip + cfg.medll_tau_step_chip * 0.5,
        cfg.medll_tau_step_chip,
    )
    separation_grid = np.arange(
        cfg.medll_separation_min_chip,
        cfg.medll_separation_max_chip + cfg.medll_separation_step_chip * 0.5,
        cfg.medll_separation_step_chip,
    )
    coarse: list[tuple[float, float, float]] = []
    for separation in separation_grid:
        costs = _batch_linear_costs(template, x, y, tau_grid, float(separation))
        coarse.extend(
            (float(cost), float(tau), float(separation))
            for tau, cost in zip(tau_grid, costs)
        )
    coarse.sort(key=lambda item: item[0])
    starts: list[tuple[float, float]] = []
    for _, tau, separation in coarse:
        if all(abs(tau - old_tau) > cfg.medll_tau_step_chip or
               abs(separation - old_sep) > cfg.medll_separation_step_chip
               for old_tau, old_sep in starts):
            starts.append((tau, separation))
        if len(starts) >= cfg.medll_refine_starts:
            break

    bounds = [
        (cfg.medll_tau_min_chip, cfg.medll_tau_max_chip),
        (cfg.medll_separation_min_chip, cfg.medll_separation_max_chip),
    ]
    best: tuple[float, np.ndarray, np.ndarray, bool] | None = None
    for start in starts:
        result = minimize(
            lambda p: _linear_fit(template, x, y, float(p[0]), float(p[1]))[0],
            np.asarray(start, dtype=float),
            method="L-BFGS-B",
            bounds=bounds,
            options={"maxiter": cfg.medll_max_nfev, "ftol": 1.0e-12},
        )
        cost, coefficient, residual = _linear_fit(
            template, x, y, float(result.x[0]), float(result.x[1])
        )
        candidate = (cost, result.x.copy(), coefficient.copy(), bool(result.success))
        if best is None or candidate[0] < best[0]:
            best = candidate
    if best is None:
        raise RuntimeError("MEDLL coarse search produced no candidate")
    cost, parameters, coefficient, success = best
    tau, separation = map(float, parameters)
    a0, a1 = coefficient[:2]
    signal_power = float(np.vdot(y, y).real)
    ratio = float(abs(a1) / max(abs(a0), 1.0e-30))
    phase = float(np.angle(a1 / a0)) if abs(a0) > 0.0 and abs(a1) > 0.0 else float("nan")
    rmse = float(np.sqrt(cost / (2.0 * len(y))))
    return {
        "tau_chip": tau,
        "separation_chip": separation,
        "amplitude_ratio": ratio,
        "relative_phase_rad": phase,
        "cost": cost,
        "fit_nrmse": float(rmse / max(np.sqrt(signal_power / (2.0 * len(y))), 1.0e-30)),
        "boundary": bool(
            np.isclose(tau, cfg.medll_tau_min_chip)
            or np.isclose(tau, cfg.medll_tau_max_chip)
            or np.isclose(separation, cfg.medll_separation_min_chip)
            or np.isclose(separation, cfg.medll_separation_max_chip)
        ),
        "optimizer_success": success,
    }


def run_baselines(
    scans: list[EpochScan],
    chain_templates: dict[int, ChainBankTemplate],
    cfg: BaselineConfig,
) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    for scan in scans:
        eml = eml_estimate(scan, cfg)
        dd = double_delta_estimate(scan, cfg)
        chain = chain_templates.get(scan.bank_window)
        if chain is None:
            raise KeyError(f"no response bank for bank_window={scan.bank_window}")
        ideal = IdealTemplate(scan.delay_chip, chain.center_chip())
        medll_ideal = medll_estimate(scan, ideal, cfg)
        medll_chain = medll_estimate(scan, chain, cfg)
        row: dict[str, object] = {
            "window": scan.window,
            "epoch": scan.epoch,
            "bank_window": scan.bank_window,
            "prompt_center_chip": ideal.center_chip,
            "eml_tau_chip": eml["tau_chip"],
            "eml_boundary": eml["boundary"],
            "eml_discriminator_residual": eml["discriminator_residual"],
            "dd_tau_chip": dd["tau_chip"],
            "dd_boundary": dd["boundary"],
            "dd_discriminator_residual": dd["discriminator_residual"],
            "medll_ideal_tau_chip": medll_ideal["tau_chip"],
            "medll_ideal_separation_chip": medll_ideal["separation_chip"],
            "medll_ideal_amplitude_ratio": medll_ideal["amplitude_ratio"],
            "medll_ideal_phase_rad": medll_ideal["relative_phase_rad"],
            "medll_ideal_cost": medll_ideal["cost"],
            "medll_ideal_fit_nrmse": medll_ideal["fit_nrmse"],
            "medll_ideal_boundary": medll_ideal["boundary"],
            "medll_ideal_optimizer_success": medll_ideal["optimizer_success"],
            "medll_chain_tau_chip": medll_chain["tau_chip"],
            "medll_chain_separation_chip": medll_chain["separation_chip"],
            "medll_chain_amplitude_ratio": medll_chain["amplitude_ratio"],
            "medll_chain_phase_rad": medll_chain["relative_phase_rad"],
            "medll_chain_cost": medll_chain["cost"],
            "medll_chain_fit_nrmse": medll_chain["fit_nrmse"],
            "medll_chain_boundary": medll_chain["boundary"],
            "medll_chain_optimizer_success": medll_chain["optimizer_success"],
        }
        output.append(row)
    return output


def _load_config(path: Path | None) -> BaselineConfig:
    if path is None:
        return BaselineConfig()
    payload = json.loads(path.read_text(encoding="utf-8"))
    allowed = set(BaselineConfig.__dataclass_fields__)
    unknown = sorted(set(payload) - allowed)
    if unknown:
        raise ValueError(f"unknown baseline config keys: {unknown}")
    return BaselineConfig(**payload)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scan-csv", type=Path, required=True)
    parser.add_argument("--chain-bank-dir", type=Path, required=True)
    parser.add_argument("--output-csv", type=Path, required=True)
    parser.add_argument("--config", type=Path)
    parser.add_argument("--expected-taps", type=int, default=401)
    parser.add_argument(
        "--bank-window", type=int,
        help="bank ID when the scanner CSV has no bank_window column",
    )
    args = parser.parse_args()
    cfg = _load_config(args.config)
    scans = read_epoch_scans(args.scan_csv, args.expected_taps, args.bank_window)
    bank_windows = sorted({scan.bank_window for scan in scans})
    templates = {
        window: ChainBankTemplate(args.chain_bank_dir / f"bank_{window}.npz")
        for window in bank_windows
    }
    rows = run_baselines(scans, templates, cfg)
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_csv.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps({
        "scan_csv": str(args.scan_csv),
        "output_csv": str(args.output_csv),
        "epochs": len(rows),
        "bank_windows": bank_windows,
        "configuration": cfg.__dict__,
        "truth_accessed": False,
    }, indent=2))


if __name__ == "__main__":
    main()
