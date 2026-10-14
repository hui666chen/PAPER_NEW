#!/usr/bin/env python3
"""Noise-whitened short-delay detectability audit using retained responses only."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq, minimize_scalar
from scipy.stats import chi2, ncx2


PAPER = Path("/home/zoxoai/Documents/ChatGPT/paper")
REPO = PAPER / "PAPER_NEW"
DATA = PAPER / "gps_l1_prn09_short_delay_pilot_20261013"
CAL = DATA / "calibration_v2_after_acquisition_amendment"
OUT = REPO / "results/short_delay_detectability_20261013"
BANKS = {"d005_a025": 39, "d005_a050": 53, "d010_a025": 67, "d010_a050": 81}
PILOT_PARAMETERS = {
    "d005_a025": (0.05, 0.25), "d005_a050": (0.05, 0.50),
    "d010_a025": (0.10, 0.25), "d010_a050": (0.10, 0.50),
}
CAL_WINDOWS = {
    "d005_a025": tuple(range(36, 44)),
    "d005_a050": tuple(range(50, 58)),
    "d010_a025": tuple(range(64, 72)),
    "d010_a050": tuple(range(78, 86)),
}
DELAYS = (0.02, 0.05, 0.10, 0.20)
AMPLITUDES = (0.10, 0.25, 0.50)
UNIFORM_PHASES_DEG = tuple(float(x) for x in np.arange(0.0, 360.0, 22.5))
PILOT_PHASE_DEG = 8.83
N_EPOCHS_AVERAGED = 100
EPOCH_SECONDS = 0.001
OBSERVATION_SECONDS = N_EPOCHS_AVERAGED * EPOCH_SECONDS
TAP_SPAN = 2.0
DELAY_GRID = np.arange(-0.40, 0.4000001, 0.001)
HAC_LAG = 20
SHRINKAGE_FLOOR = 1.0e-8
CHIP_METRES = 299_792_458.0 / 1_023_000.0

sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER / "gps_l1_prn30_formal_pvt"))
sys.path.insert(0, str(REPO / "scripts"))
from receiver_chain_response_control import ResponseSpec, load_response_model  # noqa: E402


SPEC = ResponseSpec(prn=9)
TAP = SPEC.tap_grid.astype(complex)
ONE = np.ones(SPEC.tap_count, dtype=complex)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_epoch_matrix(path: Path, wanted_windows: tuple[int, ...]) -> tuple[dict[int, np.ndarray], dict[int, dict[int, tuple[float, int]]]]:
    wanted = set(wanted_windows)
    values: dict[int, dict[int, dict[int, complex]]] = {w: {} for w in wanted_windows}
    meta: dict[int, dict[int, tuple[float, int]]] = {w: {} for w in wanted_windows}
    with path.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            window = int(row["window"])
            if window not in wanted:
                continue
            epoch, tap = int(row["epoch"]), int(row["tap"])
            values[window].setdefault(epoch, {})[tap] = complex(
                float(row["mixed_i"]), float(row["mixed_q"])
            )
            meta[window][epoch] = (float(row["time_ms"]), int(row["input_start_sample"]))
    matrices = {}
    for window in wanted_windows:
        epochs = sorted(values[window])
        if len(epochs) != 100 or epochs != list(range(100)):
            raise ValueError(f"{path}: window {window} does not have epochs 0..99")
        matrix = []
        for epoch in epochs:
            taps = values[window][epoch]
            if sorted(taps) != list(range(SPEC.tap_count)):
                raise ValueError(f"{path}: window={window} epoch={epoch} missing/duplicate taps")
            matrix.append([taps[t] for t in range(SPEC.tap_count)])
        matrices[window] = np.asarray(matrix, dtype=complex)
    return matrices, meta


def fit_clean_epoch_residuals(model, matrices: dict[int, np.ndarray]) -> tuple[list[np.ndarray], list[dict]]:
    """Remove per-epoch complex gain, small delay, and affine tap background."""
    tau_grid = np.arange(-0.20, 0.2000001, 0.01)
    nuisance = np.column_stack((ONE, TAP))
    q, _ = np.linalg.qr(nuisance, mode="reduced")
    response = model.exact(-tau_grid).T
    response_perp = response - q @ (q.conj().T @ response)
    response_norm = np.sum(np.abs(response_perp) ** 2, axis=0)
    residual_sequences = []
    window_stats = []
    for window, observations in matrices.items():
        y_perp = observations - (observations @ q) @ q.conj().T
        energy = np.sum(np.abs(y_perp) ** 2, axis=1)
        correlations = y_perp @ response_perp.conj()
        rss = energy[:, None] - np.abs(correlations) ** 2 / np.maximum(response_norm, 1e-30)
        winners = np.argmin(rss, axis=1)
        residual = np.empty_like(y_perp)
        for epoch, winner in enumerate(winners):
            p = response_perp[:, winner]
            gain = correlations[epoch, winner] / max(response_norm[winner], 1e-30)
            residual[epoch] = y_perp[epoch] - gain * p
        # Remove a deterministic within-window residual mean, not the epoch-to-epoch noise.
        residual -= residual.mean(axis=0, keepdims=True)
        residual_sequences.append(residual)
        window_stats.append({
            "window": int(window),
            "tau_fit_median_chip": float(np.median(tau_grid[winners])),
            "residual_rms_per_tap": float(np.sqrt(np.mean(np.abs(residual) ** 2))),
            "mean_signal_norm": float(np.mean(np.linalg.norm(observations, axis=1))),
        })
    return residual_sequences, window_stats


def complex_to_real_vector(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=complex)
    return np.concatenate((values.real, values.imag), axis=0)


def complex_to_real_design(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=complex)
    return np.block([[values.real, -values.imag], [values.imag, values.real]])


def complex_templates_to_real(values: np.ndarray) -> np.ndarray:
    """Map complex p×n templates to real 2p×n×2 coefficient designs."""
    values = np.asarray(values, dtype=complex)
    p, n = values.shape
    result = np.empty((2 * p, n, 2), dtype=float)
    result[:p, :, 0] = values.real
    result[p:, :, 0] = values.imag
    result[:p, :, 1] = -values.imag
    result[p:, :, 1] = values.real
    return result


def complex_jacobian_to_real(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=complex)
    return np.vstack((values.real, values.imag))


def complex_shrink_covariance(residuals: np.ndarray) -> tuple[np.ndarray, dict]:
    complex_x = np.asarray(residuals, dtype=complex)
    complex_x = complex_x - complex_x.mean(axis=0, keepdims=True)
    x = np.concatenate((complex_x.real, complex_x.imag), axis=1)
    x = x - x.mean(axis=0, keepdims=True)
    n = len(x)
    sample = (x.T @ x) / n
    target = np.diag(np.diag(sample))
    denom = float(np.linalg.norm(sample - target, "fro") ** 2)
    row_energy = np.sum(np.abs(x) ** 2, axis=1)
    phi = float((np.sum(row_energy**2) - n * np.linalg.norm(sample, "fro") ** 2) / (n**2))
    alpha = float(np.clip(max(phi, 0.0) / max(denom, 1e-30), 0.0, 1.0))
    covariance = (1.0 - alpha) * sample + alpha * target
    eigenvalues, eigenvectors = np.linalg.eigh((covariance + covariance.T) * 0.5)
    floor = max(float(np.max(eigenvalues)) * SHRINKAGE_FLOOR, 1e-20)
    clipped = np.maximum(eigenvalues, floor)
    covariance = (eigenvectors * clipped) @ eigenvectors.T
    complex_cov = (complex_x.T @ complex_x.conj()) / len(complex_x)
    pseudo = (complex_x.T @ complex_x) / len(complex_x)
    return covariance, {
        "sample_vectors": int(n), "dimension": int(x.shape[1]),
        "covariance_representation": "augmented_real_[Re(taps),Im(taps)]",
        "diagonal_shrinkage": alpha, "eigenvalue_floor": floor,
        "min_eigenvalue_after_floor": float(clipped.min()),
        "max_eigenvalue": float(clipped.max()),
        "condition_number": float(clipped.max() / clipped.min()),
        "complex_pseudo_to_covariance_frobenius_ratio": float(
            np.linalg.norm(pseudo, "fro") / max(np.linalg.norm(complex_cov, "fro"), 1e-30)
        ),
        "trace": float(np.trace(covariance).real),
    }


def temporal_diagnostics(sequences: list[np.ndarray], covariance_1ms: np.ndarray,
                         max_lag: int = HAC_LAG) -> tuple[list[dict], np.ndarray, dict]:
    lag_rows = []
    gamma_matrices = {lag: np.zeros_like(covariance_1ms) for lag in range(1, max_lag + 1)}
    gamma_traces = {lag: 0.0 for lag in range(1, max_lag + 1)}
    for x in sequences:
        x = x - x.mean(axis=0, keepdims=True)
        x = np.concatenate((x.real, x.imag), axis=1)
        for lag in range(1, max_lag + 1):
            left, right = x[:-lag], x[lag:]
            gamma = left.T @ right.conj() / len(sequences) / (100 - lag)
            gamma_matrices[lag] += gamma
            gamma_traces[lag] += float(np.trace(gamma))
    trace0 = float(np.trace(covariance_1ms).real)
    hac = covariance_1ms.astype(float).copy()
    for lag in range(1, max_lag + 1):
        rho = gamma_traces[lag] / max(trace0, 1e-30)
        bartlett = 1.0 - lag / (max_lag + 1.0)
        finite_n = 1.0 - lag / 100.0
        hac += bartlett * finite_n * (gamma_matrices[lag] + gamma_matrices[lag].conj().T)
        lag_rows.append({
            "lag_epochs": lag, "lag_seconds": lag * EPOCH_SECONDS,
            "rho_real": float(rho), "rho_imag": 0.0,
            "rho_abs": float(abs(rho)), "sequence_count": len(sequences),
        })
    hac /= 100.0
    hac = (hac + hac.T) * 0.5
    eigenvalues, eigenvectors = np.linalg.eigh(hac)
    min_before = float(eigenvalues.min())
    floor = max(float(eigenvalues.max()) * SHRINKAGE_FLOOR, 1e-20)
    clipped = np.maximum(eigenvalues, floor)
    hac = (eigenvectors * clipped) @ eigenvectors.T
    metadata = {
        "hac_lag_bandwidth": max_lag,
        "hac_eigenvalue_min_before_floor": min_before,
        "hac_eigenvalue_floor": floor,
        "hac_eigenvalue_clipped_count": int(np.count_nonzero(eigenvalues < floor)),
        "hac_trace_inflation_vs_iid": float(np.trace(hac).real / max(np.trace(covariance_1ms).real / 100.0, 1e-30)),
        "hac_condition_number": float(clipped.max() / clipped.min()),
    }
    return lag_rows, hac, metadata


def whiten_matrix(covariance: np.ndarray) -> tuple[np.ndarray, dict]:
    covariance = np.asarray(covariance.real, dtype=float)
    covariance = (covariance + covariance.T) * 0.5
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    floor = max(float(eigenvalues.max()) * SHRINKAGE_FLOOR, 1e-20)
    eigenvalues = np.maximum(eigenvalues, floor)
    whitening = (eigenvectors * (1.0 / np.sqrt(eigenvalues))) @ eigenvectors.T
    return whitening, {
        "whitening_floor": floor,
        "whitening_condition_number": float(eigenvalues.max() / eigenvalues.min()),
    }


def threshold_for_reference_df(df: int, alpha: float = 1e-3, power: float = 0.90) -> dict:
    critical = float(chi2.ppf(1.0 - alpha, df))
    nc = float(brentq(lambda value: ncx2.sf(critical, df, value) - power, 1e-12, 1e4))
    return {
        "reference_real_df": df, "false_alarm_probability": alpha,
        "target_detection_probability": power, "central_threshold": critical,
        "required_real_noncentrality": nc,
    }


def model_mean(model, delay_chip: float, amplitude: float, phase_deg: float) -> np.ndarray:
    rho = amplitude * np.exp(1j * np.deg2rad(phase_deg))
    response = model.exact([0.0, delay_chip])
    return response[0] + rho * response[1]


def make_profile_context(model, whitening: np.ndarray) -> dict:
    nuisance = complex_to_real_design(np.column_stack((ONE, TAP)))
    nuisance_w = whitening @ nuisance
    q, _ = np.linalg.qr(nuisance_w, mode="reduced")
    templates = model.exact(-DELAY_GRID).T
    template_designs = complex_templates_to_real(templates)
    templates_w = whitening @ template_designs.reshape(template_designs.shape[0], -1)
    templates_w = templates_w.reshape(template_designs.shape)
    flat_templates = templates_w.reshape(templates_w.shape[0], -1)
    projection = (q @ (q.T @ flat_templates)).reshape(templates_w.shape)
    templates_perp = templates_w - projection
    p0, p1 = templates_perp[:, :, 0], templates_perp[:, :, 1]
    g00 = np.sum(p0 * p0, axis=0)
    g11 = np.sum(p1 * p1, axis=0)
    g01 = np.sum(p0 * p1, axis=0)
    return {
        "model": model, "whitening": whitening, "q": q,
        "templates_perp": templates_perp,
        "gram": (g00, g01, g11),
        "interpolator": CubicSpline(DELAY_GRID, templates_perp.transpose(1, 0, 2), axis=0),
    }


def profile_single_path_lambda(context: dict, mean: np.ndarray) -> tuple[float, float]:
    model = context["model"]
    whitening = context["whitening"]
    q = context["q"]
    templates_perp = context["templates_perp"]
    g00, g01, g11 = context["gram"]
    interpolator = context["interpolator"]
    y_w = whitening @ complex_to_real_vector(mean)
    y_perp = y_w - q @ (q.T @ y_w)
    p0, p1 = templates_perp[:, :, 0], templates_perp[:, :, 1]
    b0, b1 = p0.T @ y_perp, p1.T @ y_perp
    determinant = np.maximum(g00 * g11 - g01 * g01, 1e-30)
    explained = (g11 * b0 * b0 + g00 * b1 * b1 - 2.0 * g01 * b0 * b1) / determinant
    scores = float(y_perp @ y_perp) - explained
    best = int(np.argmin(scores))
    lo = float(DELAY_GRID[max(0, best - 1)])
    hi = float(DELAY_GRID[min(len(DELAY_GRID) - 1, best + 1)])
    def score(tau: float) -> float:
        p = interpolator(tau)
        gram = p.T @ p
        rhs = p.T @ y_perp
        beta = np.linalg.pinv(gram, rcond=1e-12) @ rhs
        return float(y_perp @ y_perp - rhs @ beta)

    result = minimize_scalar(score, bounds=(lo, hi), method="bounded", options={"xatol": 1e-9})
    tau_hat = float(result.x if result.success else DELAY_GRID[best])
    exact_complex = model.exact([-tau_hat])[0][:, None]
    exact = whitening @ complex_to_real_design(exact_complex)
    exact_perp = exact - q @ (q.T @ exact)
    beta = np.linalg.lstsq(exact_perp, y_perp, rcond=1e-12)[0]
    residual = y_perp - exact_perp @ beta
    return float(residual @ residual), tau_hat


def local_delay_fim_width(model, delay: float, amplitude: float, phase_deg: float,
                          whitening: np.ndarray) -> tuple[float, int, float]:
    rho = amplitude * np.exp(1j * np.deg2rad(phase_deg))
    h = 1e-4

    def path(tau: float, separation: float, ratio: complex) -> np.ndarray:
        r = model.exact([-tau, separation - tau])
        return r[0] + ratio * r[1]

    base = path(0.0, delay, rho)
    d_tau = (path(h, delay, rho) - path(-h, delay, rho)) / (2.0 * h)
    d_sep = (path(0.0, delay + h, rho) - path(0.0, delay - h, rho)) / (2.0 * h)
    refl = model.exact([delay])[0]
    jac = np.column_stack((
        d_tau, d_sep, refl, 1j * refl,
        base, 1j * base, ONE, 1j * ONE, TAP, 1j * TAP,
    ))
    jw = whitening @ complex_jacobian_to_real(jac)
    column_norms = np.linalg.norm(jw, axis=0)
    normalized = jw / np.maximum(column_norms, 1e-30)[None, :]
    singular = np.linalg.svd(normalized, compute_uv=False)
    rank = int(np.count_nonzero(singular > max(float(singular.max()), 1e-30) * 1e-10))
    condition = float(singular.max() / max(singular[rank - 1], 1e-30)) if rank else float("inf")
    nuisance = jw[:, 1:] / np.maximum(column_norms[None, 1:], 1e-30)
    nuisance_u, nuisance_s, _ = np.linalg.svd(nuisance, full_matrices=False)
    nuisance_rank = int(np.count_nonzero(nuisance_s > max(float(nuisance_s.max()), 1e-30) * 1e-10))
    qn = nuisance_u[:, :nuisance_rank]
    tau_residual = jw[:, 0] - qn @ (qn.T @ jw[:, 0])
    efficient = float(tau_residual @ tau_residual)
    sigma_chip = float(1.0 / np.sqrt(efficient)) if efficient > 1e-20 else float("inf")
    width_m = 2.0 * 1.959963984540054 * sigma_chip * CHIP_METRES
    if rank < jac.shape[1]:
        width_m = float("inf")
    return width_m, rank, condition


def write_csv(path: Path, rows: list[dict]) -> None:
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    protocol = OUT / "FROZEN_SHORT_DELAY_DETECTABILITY_20261013.md"
    if not protocol.is_file():
        raise FileNotFoundError(protocol)
    thresholds = {df: threshold_for_reference_df(df) for df in (1, 3)}
    all_rows = []
    covariance_summary = {}
    lag_records = []
    input_hashes = {str(protocol): sha256(protocol)}
    for cell, bank_window in BANKS.items():
        source = CAL / "open_loop" / cell / "prn09_epoch_responses.csv"
        bank_path = CAL / "response_banks" / f"bank_{bank_window}.npz"
        model = load_response_model(bank_path, expected_spec=SPEC)
        matrices, metadata = read_epoch_matrix(source, CAL_WINDOWS[cell])
        residual_sequences, window_stats = fit_clean_epoch_residuals(model, matrices)
        residuals = np.concatenate(residual_sequences, axis=0)
        covariance_1ms, cov_meta = complex_shrink_covariance(residuals)
        lag_rows, covariance_hac, hac_meta = temporal_diagnostics(residual_sequences, covariance_1ms)
        for item in lag_rows:
            lag_records.append({"cell": cell, **item})
        covariance_iid = covariance_1ms / N_EPOCHS_AVERAGED
        lag_abs_max = max(item["rho_abs"] for item in lag_rows)
        whitening_iid, whiten_meta_iid = whiten_matrix(covariance_iid)
        whitening_hac, whiten_meta_hac = whiten_matrix(covariance_hac)
        profile_iid = make_profile_context(model, whitening_iid)
        profile_hac = make_profile_context(model, whitening_hac)
        covariance_summary[cell] = {
            **cov_meta,
            "calibration_windows": list(CAL_WINDOWS[cell]),
            "epochs_per_window": 100,
            "epoch_duration_seconds": EPOCH_SECONDS,
            "averaged_span_seconds": OBSERVATION_SECONDS,
            "window_start_spacing_seconds": 1.0,
            "lag1_to20_abs_max": lag_abs_max,
            **hac_meta,
            "iid_whitening": whiten_meta_iid,
            "hac_whitening": whiten_meta_hac,
            "per_window_fit": window_stats,
        }
        input_hashes[str(source)] = sha256(source)
        input_hashes[str(bank_path)] = sha256(bank_path)
        phases = tuple(sorted(set(UNIFORM_PHASES_DEG + (PILOT_PHASE_DEG,))))
        for delay in DELAYS:
            for amp in AMPLITUDES:
                for phase in phases:
                    mean = model_mean(model, delay, amp, phase)
                    for covariance_kind, whitening, profile_context in (
                        ("iid_div100", whitening_iid, profile_iid),
                        ("hac", whitening_hac, profile_hac),
                    ):
                        lambda_real, tau_hat = profile_single_path_lambda(profile_context, mean)
                        width95_m, fim_rank, fim_condition = local_delay_fim_width(
                            model, delay, amp, phase, whitening
                        )
                        row = {
                            "cell": cell,
                            "bank_window": bank_window,
                            "delay_chip": delay,
                            "amplitude_ratio": amp,
                            "phase_deg": phase,
                            "phase_kind": "pilot_actual" if phase == PILOT_PHASE_DEG else "uniform_grid",
                            "covariance_model": covariance_kind,
                            "observed_integration_seconds": OBSERVATION_SECONDS,
                            "single_path_tau_hat_chip": tau_hat,
                            "lambda_augmented_real_rT_Cinv_r": lambda_real,
                            "lambda_real_noncentrality": lambda_real,
                            "lambda_real_per_second_extrapolation": lambda_real / OBSERVATION_SECONDS,
                            "time_to_df1_threshold_seconds": OBSERVATION_SECONDS * thresholds[1]["required_real_noncentrality"] / max(lambda_real, 1e-300),
                            "time_to_df3_threshold_seconds": OBSERVATION_SECONDS * thresholds[3]["required_real_noncentrality"] / max(lambda_real, 1e-300),
                            "passes_df1_reference": lambda_real >= thresholds[1]["required_real_noncentrality"],
                            "passes_df3_reference": lambda_real >= thresholds[3]["required_real_noncentrality"],
                            "local_fim_95_direct_delay_width_m": width95_m,
                            "fim_effective_rank": fim_rank,
                            "fim_condition_number": fim_condition,
                            "eligible_for_existing_data_refinement": bool(
                                phase == PILOT_PHASE_DEG
                                and covariance_kind == "hac"
                                and delay == PILOT_PARAMETERS[cell][0]
                                and amp == PILOT_PARAMETERS[cell][1]
                                and lambda_real >= thresholds[3]["required_real_noncentrality"]
                                and width95_m <= 3.0
                            ),
                        }
                        all_rows.append(row)
    write_csv(OUT / "detectability_grid.csv", all_rows)
    write_csv(OUT / "temporal_lag_diagnostics.csv", lag_records)
    summary = {
        "status": "complete",
        "analysis_class": "post_hoc_theoretical_diagnostic_existing_data_only",
        "iq_generated": False,
        "pvt_run": False,
        "timing": {
            "stored_epoch_timestamp_step_ms": 1.0,
            "stored_epoch_sample_step": 25000,
            "sample_rate_hz": 25_000_000,
            "correlator_epoch_seconds": EPOCH_SECONDS,
            "100_epoch_average_span_seconds": OBSERVATION_SECONDS,
            "one_second_start_spacing": True,
            "warning": "The retained 100-epoch average spans 100 ms, not a full 1 s; 1-s time-to-detection is an extrapolation.",
        },
        "threshold_references": thresholds,
        "covariance_by_cell": covariance_summary,
        "row_count": len(all_rows),
        "uniform_phase_points": len(UNIFORM_PHASES_DEG),
        "pilot_phase_points": 1,
        "eligible_rows": sum(bool(row["eligible_for_existing_data_refinement"]) for row in all_rows),
        "eligibility_rule": "actual pilot phase 8.83 deg; HAC covariance; df=3 reference threshold; local 95% FIM width <=3m",
        "input_sha256": input_hashes,
    }
    (OUT / "SUMMARY.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    pilot_rows = [
        row for row in all_rows
        if row["phase_kind"] == "pilot_actual"
        and row["covariance_model"] == "hac"
        and row["delay_chip"] == PILOT_PARAMETERS[row["cell"]][0]
        and row["amplitude_ratio"] == PILOT_PARAMETERS[row["cell"]][1]
    ]
    report = [
        "# Short-delay detectability diagnostic",
        "",
        "Existing calibration and validation response data only; no IQ was generated and no PVT was run.",
        "",
        "## Timing audit",
        "",
        "The stored per-epoch timestamps advance by 1 ms and sample starts by 25,000 samples at 25 MHz. Thus each 100-epoch averaged 401-tap profile spans 100 ms. The nominal profile starts are one second apart; this dataset does not contain a full 1-s coherent average. Any 1-s time-to-detection below is a linear, same-phase extrapolation from 100 ms.",
        "",
        "## Reference thresholds",
        "",
        f"- df=1: required real noncentrality {thresholds[1]['required_real_noncentrality']:.3f}.",
        f"- df=3: required real noncentrality {thresholds[3]['required_real_noncentrality']:.3f}.",
        "- These are chi-square reference thresholds for PFA=1e-3 and power=0.90, not an exact ROC calibration for the non-regular unknown-delay multipath GLRT.",
        "",
        "## Pilot cells at the configured 8.83° phase (HAC covariance)",
        "",
        "| Cell | λ real | Time to df=3 threshold (s) | Local 95% delay width (m) | Eligible for dense reanalysis |",
        "|---|---:|---:|---:|---:|",
    ]
    for row in pilot_rows:
        report.append(
            f"| {row['cell']} | {row['lambda_real_noncentrality']:.4g} | {row['time_to_df3_threshold_seconds']:.4g} | {row['local_fim_95_direct_delay_width_m']:.4g} | {row['eligible_for_existing_data_refinement']} |"
        )
    report += [
        "",
        f"Eligible pilot rows under the frozen rule: {summary['eligible_rows']}.",
        "",
        "Full phase sweep, covariance diagnostics, per-window calibration fit diagnostics, and all λ/FIM values are in `detectability_grid.csv`, `temporal_lag_diagnostics.csv`, and `SUMMARY.json`.",
    ]
    (OUT / "REPORT.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": summary["status"], "rows": len(all_rows),
        "eligible_rows": summary["eligible_rows"],
        "pilot_hac": [{
            "cell": row["cell"], "lambda_real": row["lambda_real_noncentrality"],
            "width95_m": row["local_fim_95_direct_delay_width_m"],
            "time_to_df3_s": row["time_to_df3_threshold_seconds"],
        } for row in pilot_rows],
    }, indent=2), flush=True)


if __name__ == "__main__":
    main()
