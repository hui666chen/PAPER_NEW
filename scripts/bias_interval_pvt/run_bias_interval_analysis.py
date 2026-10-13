#!/usr/bin/env python3
"""Post-hoc bias-interval reanalysis of retained PRN09 correlator data.

No IQ is generated. Every validation grid residual is retained in compressed
NPZ form. The short-delay existence gate is recomputed observation-only using
the already frozen threshold and 3-of-5 rule.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

import numpy as np


WORKSPACE = Path("/home/zoxoai/Documents/ChatGPT/paper")
REPO = WORKSPACE / "PAPER_NEW"
DATA_LONG = WORKSPACE / "gps_l1_prn09_comparable_v3"
DATA_SHORT = WORKSPACE / "gps_l1_prn09_short_delay_pilot_20261013"
PRN30 = WORKSPACE / "gps_l1_prn30_formal_pvt"
RESULTS = REPO / "results/bias_interval_pvt_20261013"
SURFACES = RESULTS / "residual_surfaces"
CHIP_METRES = 299_792_458.0 / 1_023_000.0
BOOTSTRAP_SEED = 2026101316
BOOTSTRAPS_PER_CLEAN_WINDOW = 100
CALIBRATION_QUANTILE = 0.90
VALIDATION_SEEDS_LONG = (2026101110, 2026101111, 2026101112, 2026101113, 2026101114)
VALIDATION_SEEDS_SHORT = (2026101311, 2026101312, 2026101313, 2026101314, 2026101320)
LONG_SCENES = {
    "clean": {"windows": tuple(range(61, 77)), "bank_window": 68},
    "amp_0p05": {"windows": tuple(range(100, 116)), "bank_window": 108},
    "amp_0p10": {"windows": tuple(range(121, 137)), "bank_window": 129},
}
SHORT_CELLS = {
    "d005_a025": {"windows": tuple(range(36, 44)), "bank_window": 39, "multipath": True},
    "d005_a050": {"windows": tuple(range(50, 58)), "bank_window": 53, "multipath": True},
    "d010_a025": {"windows": tuple(range(64, 72)), "bank_window": 67, "multipath": True},
    "d010_a050": {"windows": tuple(range(78, 86)), "bank_window": 81, "multipath": True},
}

sys.path.insert(0, str(WORKSPACE))
sys.path.insert(0, str(PRN30))
sys.path.insert(0, str(WORKSPACE / "gps_l1_prn09_shared_ratio_v2"))
sys.path.insert(0, str(REPO / "scripts/prn09_comparable_v3"))

import shared_ratio_existence as shared  # noqa: E402
from monte_carlo_resolution_baseline import build_dictionary  # noqa: E402
from receiver_chain_response_control import ResponseSpec, load_response_model  # noqa: E402
from estimate_formal_windows import load_windows  # noqa: E402


SPEC = ResponseSpec(prn=9)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def higher_quantile(values: np.ndarray, probability: float) -> float:
    """Order-statistic quantile using the higher-neighbor convention."""
    ordered = np.sort(np.asarray(values, dtype=float))
    if not ordered.size:
        raise ValueError("cannot compute a quantile from an empty sample")
    index = int(np.ceil(float(probability) * (len(ordered) - 1)))
    return float(ordered[index])


def read_average(path: Path) -> dict[int, dict[str, np.ndarray]]:
    grouped: dict[int, dict[str, list]] = {}
    with path.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            window = int(row["window"])
            item = grouped.setdefault(window, {"delay": [], "mixed": []})
            item["delay"].append(float(row["delay_chip"]))
            item["mixed"].append(complex(float(row["mixed_i"]), float(row["mixed_q"])))
    result = {}
    for window, item in grouped.items():
        delay = np.asarray(item["delay"], dtype=float)
        y = np.asarray(item["mixed"], dtype=complex)
        if len(delay) != 401 or len(y) != 401 or np.any(np.diff(delay) <= 0):
            raise ValueError(f"{path}: window {window} does not contain 401 ordered taps")
        result[window] = {"delay": delay, "y": y}
    return result


def read_epoch_matrix(path: Path, wanted_windows: tuple[int, ...]) -> dict[int, np.ndarray]:
    rows: dict[tuple[int, int], list[complex]] = {}
    wanted = set(wanted_windows)
    with path.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            window = int(row["window"])
            if window not in wanted:
                continue
            epoch = int(row["epoch"])
            rows.setdefault((window, epoch), []).append(
                complex(float(row["mixed_i"]), float(row["mixed_q"]))
            )
    output: dict[int, np.ndarray] = {}
    for window in wanted_windows:
        epochs = []
        for (candidate_window, epoch), values in sorted(rows.items()):
            if candidate_window == window:
                if len(values) != 401:
                    raise ValueError(f"{path}: window={window}, epoch={epoch}: {len(values)} taps")
                epochs.append(np.asarray(values, dtype=complex))
        if not epochs:
            raise ValueError(f"{path}: no epoch rows for window {window}")
        output[window] = np.stack(epochs)
    return output


class GridProfiler:
    """Fast exact variable-projection RSS over the frozen shared-ratio grid."""

    def __init__(self, bank_path: Path):
        self.bank_path = bank_path
        model = load_response_model(bank_path, expected_spec=SPEC)
        delay_grid, dictionary, _, _, _ = build_dictionary(model, SPEC)
        self.dictionary = dictionary
        self.delay_grid = delay_grid
        self.rho_grid = shared.make_rho_grid()
        self.candidate_rows, self.candidate_tau, self.candidate_sep = shared.build_candidates(
            dictionary, delay_grid
        )
        self.direct = self.candidate_rows[:, 0]
        self.reflected = self.candidate_rows[:, 1]
        gram_all = dictionary.conj().T @ dictionary
        constant_index = dictionary.shape[1] - 2
        linear_index = dictionary.shape[1] - 1
        self.inverse_grams = np.empty(
            (len(self.rho_grid), len(self.direct), 3, 3), dtype=complex
        )
        for ir, rho in enumerate(self.rho_grid):
            d, r = self.direct, self.reflected
            gram = np.zeros((len(self.direct), 3, 3), dtype=complex)
            gram[:, 0, 0] = (
                gram_all[d, d]
                + rho * gram_all[d, r]
                + np.conj(rho) * gram_all[r, d]
                + (abs(rho) ** 2) * gram_all[r, r]
            )
            gram[:, 0, 1] = gram_all[d, constant_index] + np.conj(rho) * gram_all[r, constant_index]
            gram[:, 0, 2] = gram_all[d, linear_index] + np.conj(rho) * gram_all[r, linear_index]
            gram[:, 1, 0] = np.conj(gram[:, 0, 1])
            gram[:, 2, 0] = np.conj(gram[:, 0, 2])
            gram[:, 1, 1] = gram_all[constant_index, constant_index]
            gram[:, 1, 2] = gram_all[constant_index, linear_index]
            gram[:, 2, 1] = np.conj(gram[:, 1, 2])
            gram[:, 2, 2] = gram_all[linear_index, linear_index]
            # Preserve the exact numerical rank convention of the frozen
            # shared-ratio implementation so every saved surface is comparable.
            self.inverse_grams[ir] = np.linalg.pinv(gram, rcond=1e-12)
        self.unique_tau = np.unique(self.candidate_tau)
        self.tau_groups = [np.flatnonzero(np.isclose(self.candidate_tau, tau)) for tau in self.unique_tau]
        self.grid_shape = (len(self.rho_grid), len(self.candidate_tau))

    def costs(self, y: np.ndarray) -> np.ndarray:
        y = np.asarray(y, dtype=complex)
        if y.shape != (401,):
            raise ValueError(f"expected one 401-tap vector, got {y.shape}")
        corr = y @ self.dictionary.conj()
        energy = float(np.vdot(y, y).real)
        constant_index = self.dictionary.shape[1] - 2
        linear_index = self.dictionary.shape[1] - 1
        out = np.empty(self.grid_shape, dtype=np.float64)
        d, r = self.direct, self.reflected
        for ir, rho in enumerate(self.rho_grid):
            rhs = np.empty((len(d), 3), dtype=complex)
            rhs[:, 0] = corr[d] + np.conj(rho) * corr[r]
            rhs[:, 1] = corr[constant_index]
            rhs[:, 2] = corr[linear_index]
            explained = np.einsum(
                "ni,nij,nj->n", rhs.conj(), self.inverse_grams[ir], rhs,
                optimize=True,
            ).real
            residual = energy - explained
            out[ir] = np.maximum(residual, 0.0)
        return out

    def tau_min_costs(self, costs: np.ndarray) -> np.ndarray:
        return np.asarray([float(np.min(costs[:, group])) for group in self.tau_groups])

    def required_kappa_for_zero(self, y: np.ndarray, tau_zero: float) -> float:
        costs = self.costs(y)
        tau_cost = self.tau_min_costs(costs)
        minimum = float(np.min(tau_cost))
        exact = np.flatnonzero(np.isclose(self.unique_tau, tau_zero, atol=1e-10, rtol=0.0))
        if exact.size:
            required = float(tau_cost[int(exact[0])])
        else:
            left = np.flatnonzero(self.unique_tau < tau_zero)
            right = np.flatnonzero(self.unique_tau > tau_zero)
            if not left.size or not right.size:
                return float("inf")
            required = max(float(np.min(tau_cost[left])), float(np.min(tau_cost[right])))
        return max(0.0, required / max(minimum, 1e-20) - 1.0)

    def interval(self, y: np.ndarray, tau_zero: float, kappa: float, gate_active: bool) -> dict:
        costs = self.costs(y)
        minimum = float(np.min(costs))
        near = costs <= minimum * (1.0 + float(kappa)) + 1e-14 * max(minimum, 1.0)
        if gate_active:
            selected_tau = self.candidate_tau[np.any(near, axis=0)]
            if not selected_tau.size:
                raise RuntimeError("empty retained candidate set")
            corrections = -(selected_tau - tau_zero) * CHIP_METRES
            lower, upper = float(np.min(corrections)), float(np.max(corrections))
            tau_lower, tau_upper = float(np.min(selected_tau)), float(np.max(selected_tau))
        else:
            lower = upper = tau_lower = tau_upper = 0.0
        return {
            "b_min_m": lower,
            "b_max_m": upper,
            "width_m": upper - lower,
            "tau_min_chip": tau_lower,
            "tau_max_chip": tau_upper,
            "rss_min": minimum,
            "retained_candidate_count": int(np.count_nonzero(near)) if gate_active else 0,
            "grid_candidate_count": int(costs.size),
        }, costs


def load_gate_freeze() -> tuple[float, int, int]:
    threshold_path = WORKSPACE / "gps_l1_prn30_shared_ratio_v2/null_calibration/threshold_freeze.json"
    rule_path = WORKSPACE / "gps_l1_prn30_shared_ratio_v2/PERSISTENCE_RULE_FREEZE.json"
    threshold_payload = json.loads(threshold_path.read_text(encoding="utf-8"))
    rule_payload = json.loads(rule_path.read_text(encoding="utf-8"))
    threshold = float(threshold_payload["T_null"])
    if "rule" in rule_payload:
        rule = rule_payload["rule"]
    elif "selection" in rule_payload and "selected" in rule_payload["selection"]:
        rule = rule_payload["selection"]["selected"]
    else:
        rule = rule_payload
    return threshold, int(rule["k"]), int(rule["m"])


def load_long_gate(seed: int) -> dict[tuple[str, int], bool]:
    path = DATA_LONG / f"confirmation_observation/seed_{seed}/shared_ratio/gated_records.json"
    payload = json.loads(path.read_text(encoding="utf-8"))
    return {(row["scene"], int(row["window"])): bool(row["gate_active"]) for row in payload}


def compute_short_gate(seed: int, threshold: float, k: int, m: int, profilers: dict[int, GridProfiler]) -> dict:
    output = {}
    evidence_by_cell = {}
    for cell, spec in SHORT_CELLS.items():
        run_dir = DATA_SHORT / f"seed_{seed}/open_loop/{cell}"
        epoch_path = run_dir / "prn09_epoch_responses.csv"
        windows = spec["windows"]
        observations = load_windows(epoch_path, windows)
        profiler = profilers[spec["bank_window"]]
        model = load_response_model(
            DATA_SHORT / f"calibration_v2_after_acquisition_amendment/response_banks/bank_{spec['bank_window']}.npz",
            expected_spec=SPEC,
        )
        delay_grid, dictionary, _, _, _ = build_dictionary(model, SPEC)
        rho_values = shared.make_rho_grid()
        evidence = []
        raw_records = []
        for window in windows:
            values, weights, _ = observations[window]
            record = shared.analyze_window(values, weights, dictionary, delay_grid, rho_values)
            exists = bool(record["cv_gain_min"] > threshold and record["cv_no_boundary"])
            evidence.append(exists)
            raw_records.append({
                "cell": cell,
                "window": int(window),
                "cv_gain_min": float(record["cv_gain_min"]),
                "cv_no_boundary": bool(record["cv_no_boundary"]),
                "single_window_exists": exists,
            })
        active = []
        for index in range(len(windows)):
            active.append(bool(index >= m - 1 and sum(evidence[index-m+1:index+1]) >= k))
        for row, is_active in zip(raw_records, active):
            row["gate_active"] = is_active
            output[(cell, int(row["window"]))] = is_active
        evidence_by_cell[cell] = raw_records
    gate_path = RESULTS / f"short_gate_seed_{seed}.json"
    gate_path.write_text(json.dumps({
        "seed": seed,
        "truth_accessed": False,
        "threshold": threshold,
        "rule": {"k": k, "m": m},
        "records": evidence_by_cell,
    }, indent=2) + "\n", encoding="utf-8")
    return output


def one_path_tau_zero(average_path: Path, window: int, bank_path: Path) -> float:
    averages = read_average(average_path)
    model = load_response_model(bank_path, expected_spec=SPEC)
    delay_grid, dictionary, _, _, _ = build_dictionary(model, SPEC)
    y = averages[window]["y"][None, :]
    fit = shared.one_path_profile(y, np.ones(1), dictionary, delay_grid)
    return float(fit["tau_chip"])


def calibrate_kappa(items: list[dict], tau_zero: dict[str, float], profiler_by_key: dict[str, GridProfiler]) -> dict:
    required = []
    per_window = []
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    for item in items:
        epochs = item["epochs"]
        if len(epochs) < 2:
            raise ValueError("calibration window needs at least two retained epochs")
        profiler = profiler_by_key[item["profiler_key"]]
        center = tau_zero[item["stratum"]]
        window_required = []
        for replicate in range(BOOTSTRAPS_PER_CLEAN_WINDOW):
            indices = rng.integers(0, len(epochs), size=len(epochs))
            sample = np.mean(epochs[indices], axis=0)
            kappa = profiler.required_kappa_for_zero(sample, center)
            required.append(kappa)
            window_required.append(kappa)
        per_window.append({
            "stratum": item["stratum"],
            "window": item["window"],
            "bootstrap_replicates": BOOTSTRAPS_PER_CLEAN_WINDOW,
            "required_kappa_q50": float(np.quantile(window_required, 0.50)),
            "required_kappa_q90": float(np.quantile(window_required, 0.90)),
        })
    values = np.asarray(required, dtype=float)
    if not np.all(np.isfinite(values)):
        raise RuntimeError("non-finite calibration κ requirement")
    kappa = higher_quantile(values, CALIBRATION_QUANTILE)
    return {
        "kappa": kappa,
        "calibration_quantile": CALIBRATION_QUANTILE,
        "bootstrap_seed": BOOTSTRAP_SEED,
        "bootstrap_replicates_per_window": BOOTSTRAPS_PER_CLEAN_WINDOW,
        "calibration_draw_count": int(values.size),
        "required_kappa_quantiles": {
            str(q): float(np.quantile(values, q)) for q in (0.5, 0.9, 0.95, 0.99)
        },
        "per_window": per_window,
    }


def save_surface(path: Path, costs: np.ndarray, profiler: GridProfiler, metadata: dict) -> None:
    if path.exists():
        raise FileExistsError(f"refusing to overwrite retained residual surface: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        path,
        rss=costs,
        rho_grid=profiler.rho_grid,
        candidate_tau_chip=profiler.candidate_tau,
        candidate_separation_chip=profiler.candidate_sep,
        metadata_json=np.asarray(json.dumps(metadata, sort_keys=True)),
    )


def read_truth_short(seed: int) -> dict[tuple[str, int], float]:
    path = DATA_SHORT / f"seed_{seed}/epoch_scores.csv"
    output = {}
    with path.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            if row["method"] != "raw":
                continue
            output[(row["cell"], int(row["second"]))] = float(row["raw_error_m"])
    return output


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    fields = list(rows[0])
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--skip-gate-recompute", action="store_true")
    args = parser.parse_args()
    RESULTS.mkdir(parents=True, exist_ok=True)
    SURFACES.mkdir(parents=True, exist_ok=True)
    protocol = RESULTS / "FROZEN_REANALYSIS_RULES_20261013.md"
    if not protocol.is_file():
        raise FileNotFoundError(protocol)

    long_cal = DATA_LONG / "calibration/gpu_reference"
    long_banks = DATA_LONG / "response_banks"
    short_cal = DATA_SHORT / "calibration_v2_after_acquisition_amendment"
    short_banks = short_cal / "response_banks"

    long_profiler: dict[int, GridProfiler] = {}
    for bank_window in (68, 108, 129):
        long_profiler[bank_window] = GridProfiler(long_banks / f"bank_{bank_window}.npz")
    short_profiler: dict[int, GridProfiler] = {}
    for spec in SHORT_CELLS.values():
        bank_window = int(spec["bank_window"])
        if bank_window not in short_profiler:
            short_profiler[bank_window] = GridProfiler(short_banks / f"bank_{bank_window}.npz")

    # Independent clean zero points, fitted with the exact same one-path model.
    long_tau0 = {}
    for scene, scene_spec in LONG_SCENES.items():
        window = scene_spec["bank_window"]
        avg_path = long_cal / f"window_{window}/frozen_open_loop_prn9_responses.csv"
        bank_path = long_banks / f"bank_{window}.npz"
        long_tau0[scene] = one_path_tau_zero(avg_path, window, bank_path)
    short_tau0 = {}
    for cell, cell_spec in SHORT_CELLS.items():
        windows = cell_spec["windows"]
        avg_path = short_cal / f"open_loop/{cell}/frozen_open_loop_prn9_responses.csv"
        bank_path = short_banks / f"bank_{cell_spec['bank_window']}.npz"
        values = [one_path_tau_zero(avg_path, window, bank_path) for window in windows]
        short_tau0[cell] = float(np.median(values))

    # κ is calibrated solely from epoch responses in the independent clean seeds.
    long_cal_items = []
    long_cal_prof = {}
    for scene, spec in LONG_SCENES.items():
        window = int(spec["bank_window"])
        key = f"long_{window}"
        epochs = read_epoch_matrix(
            long_cal / f"window_{window}/prn09_epoch_responses.csv", (window,)
        )[window]
        long_cal_items.append({"stratum": scene, "window": window, "epochs": epochs, "profiler_key": key})
        long_cal_prof[key] = long_profiler[window]
    long_kappa = calibrate_kappa(long_cal_items, long_tau0, long_cal_prof)

    short_cal_items = []
    short_cal_prof = {}
    for cell, spec in SHORT_CELLS.items():
        key = f"short_{cell}"
        windows = spec["windows"]
        epoch_data = read_epoch_matrix(short_cal / f"open_loop/{cell}/prn09_epoch_responses.csv", windows)
        for window in windows:
            short_cal_items.append({"stratum": cell, "window": window, "epochs": epoch_data[window], "profiler_key": key})
        short_cal_prof[key] = short_profiler[int(spec["bank_window"])]
    short_kappa = calibrate_kappa(short_cal_items, short_tau0, short_cal_prof)

    threshold, gate_k, gate_m = load_gate_freeze()
    short_gate_by_seed = {}
    if not args.skip_gate_recompute:
        for seed in VALIDATION_SEEDS_SHORT:
            print(f"gate seed={seed}", flush=True)
            short_gate_by_seed[seed] = compute_short_gate(seed, threshold, gate_k, gate_m, short_profiler)
    else:
        raise ValueError("short-delay post-hoc gate records are required; do not skip recomputation")

    interval_rows = []
    truth_rows = []
    gate_rows = []

    for seed in VALIDATION_SEEDS_LONG:
        gates = load_long_gate(seed)
        for scene, spec in LONG_SCENES.items():
            average_path = DATA_LONG / f"confirmation_observation/seed_{seed}/open_loop/{scene}/frozen_open_loop_prn9_responses.csv"
            averages = read_average(average_path)
            profiler = long_profiler[int(spec["bank_window"])]
            for window in spec["windows"]:
                gate = bool(gates[(scene, window)])
                profile, costs = profiler.interval(averages[window]["y"], long_tau0[scene], long_kappa["kappa"], gate)
                out_path = SURFACES / "long_0p8" / f"seed_{seed}" / f"{scene}_window_{window}.npz"
                save_surface(out_path, costs, profiler, {
                    "seed": seed, "scene": scene, "window": window,
                    "bank_window": spec["bank_window"], "tau_zero_chip": long_tau0[scene],
                    "kappa": long_kappa["kappa"], "gate_active": gate,
                    "rss_definition": "unweighted complex LS after profiling common complex gain and affine tap nuisance",
                })
                interval_rows.append({
                    "family": "long_0p8", "seed": seed, "scene": scene, "cell": scene,
                    "window": window, "gate_active": gate, "tau_zero_chip": long_tau0[scene],
                    "kappa": long_kappa["kappa"], **profile, "surface_file": str(out_path),
                })
                gate_rows.append({"family": "long_0p8", "seed": seed, "scene": scene, "window": window,
                                  "injected_multipath": scene != "clean", "gate_active": gate})

    for seed in VALIDATION_SEEDS_SHORT:
        gates = short_gate_by_seed[seed]
        truth = read_truth_short(seed)
        for cell, spec in SHORT_CELLS.items():
            average_path = DATA_SHORT / f"seed_{seed}/open_loop/{cell}/frozen_open_loop_prn9_responses.csv"
            averages = read_average(average_path)
            profiler = short_profiler[int(spec["bank_window"])]
            for window in spec["windows"]:
                gate = bool(gates[(cell, window)])
                profile, costs = profiler.interval(averages[window]["y"], short_tau0[cell], short_kappa["kappa"], gate)
                out_path = SURFACES / "short_delay" / f"seed_{seed}" / f"{cell}_window_{window}.npz"
                save_surface(out_path, costs, profiler, {
                    "seed": seed, "cell": cell, "window": window,
                    "bank_window": spec["bank_window"], "tau_zero_chip": short_tau0[cell],
                    "kappa": short_kappa["kappa"], "gate_active": gate,
                    "rss_definition": "unweighted complex LS after profiling common complex gain and affine tap nuisance",
                })
                interval_rows.append({
                    "family": "short_delay", "seed": seed, "scene": "injected", "cell": cell,
                    "window": window, "gate_active": gate, "tau_zero_chip": short_tau0[cell],
                    "kappa": short_kappa["kappa"], **profile, "surface_file": str(out_path),
                })
                true_raw_error = truth[(cell, window)]
                true_correction = -true_raw_error
                truth_rows.append({
                    "family": "short_delay", "seed": seed, "cell": cell, "window": window,
                    "raw_counterfactual_error_m": true_raw_error,
                    "true_correction_m": true_correction,
                    "b_min_m": profile["b_min_m"], "b_max_m": profile["b_max_m"],
                    "covered": bool(profile["b_min_m"] - 1e-9 <= true_correction <= profile["b_max_m"] + 1e-9),
                    "gate_active": gate, "width_m": profile["width_m"],
                })
                gate_rows.append({"family": "short_delay", "seed": seed, "scene": "injected",
                                  "cell": cell, "window": window, "injected_multipath": True,
                                  "gate_active": gate})

    write_csv(RESULTS / "per_window_bias_intervals.csv", interval_rows)
    write_csv(RESULTS / "short_delay_counterfactual_coverage.csv", truth_rows)
    write_csv(RESULTS / "gate_window_records.csv", gate_rows)

    summary = {
        "status": "stage1_and_stage2_existing_data_reanalysis_complete",
        "data_generated": False,
        "analysis_class": "exploratory_post_hoc_reanalysis",
        "protocol_sha256": sha256(protocol),
        "grid": {
            "tau_chip": [float(shared.TAU[0]), float(shared.TAU[-1]), 0.02],
            "separation_chip": [float(shared.SEPARATION[0]), float(shared.SEPARATION[-1]), 0.02],
            "candidate_count": int(len(long_profiler[68].candidate_tau)),
            "rho_count": int(len(shared.make_rho_grid())),
            "rho_magnitudes": [0.0, 0.01, 0.02, 0.05, 0.10, 0.20, 0.34, 0.50],
            "phase_step_deg": 22.5,
        },
        "gate": {"threshold": threshold, "k": gate_k, "m": gate_m},
        "calibration": {"long_0p8": long_kappa, "short_delay": short_kappa,
                        "tau_zero_chip": {"long_0p8": long_tau0, "short_delay": short_tau0}},
        "long_0p8": {
            "windows": int(sum(row["family"] == "long_0p8" for row in interval_rows)),
            "counterfactual_coverage": None,
            "coverage_status": "not_scoreable: no same-seed same-second clean RINEX counterfactual retained",
        },
        "short_delay": {
            "windows": len(truth_rows),
            "coverage_rate": float(np.mean([row["covered"] for row in truth_rows])),
            "median_width_m": float(np.median([row["width_m"] for row in truth_rows])),
            "mean_width_m": float(np.mean([row["width_m"] for row in truth_rows])),
            "coverage_by_cell": {},
            "coverage_by_delay_chip": {},
        },
        "gate_rates": {},
    }
    for cell in SHORT_CELLS:
        selected = [row for row in truth_rows if row["cell"] == cell]
        summary["short_delay"]["coverage_by_cell"][cell] = {
            "coverage": float(np.mean([row["covered"] for row in selected])),
            "median_width_m": float(np.median([row["width_m"] for row in selected])),
            "gate_active_rate": float(np.mean([row["gate_active"] for row in selected])),
        }
    for delay in (0.05, 0.10):
        selected = [row for row in truth_rows if row["cell"].startswith(f"d{int(delay*1000):03d}")]
        summary["short_delay"]["coverage_by_delay_chip"][str(delay)] = {
            "coverage": float(np.mean([row["covered"] for row in selected])),
            "median_width_m": float(np.median([row["width_m"] for row in selected])),
            "gate_active_rate": float(np.mean([row["gate_active"] for row in selected])),
        }
    for family, scenes in (("long_0p8", ("clean", "amp_0p05", "amp_0p10")), ("short_delay", ("injected",))):
        for scene in scenes:
            selected = [row for row in gate_rows if row["family"] == family and row["scene"] == scene]
            if not selected:
                continue
            active = float(np.mean([row["gate_active"] for row in selected]))
            is_mp = bool(selected[0]["injected_multipath"])
            summary["gate_rates"][f"{family}/{scene}"] = {
                "windows": len(selected), "active_rate": active,
                "false_positive_rate": active if not is_mp else None,
                "false_negative_rate": 1.0 - active if is_mp else None,
            }
    short_coverage = summary["short_delay"]["coverage_rate"]
    summary["stage2_decision"] = "PASS_TO_PVT" if short_coverage >= 0.90 else "STOP_BEFORE_PVT"
    (RESULTS / "SUMMARY.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    lines = [
        "# Bias-interval reanalysis results",
        "",
        "Classification: exploratory/post-hoc; no IQ was generated.",
        "",
        f"- Grid: {summary['grid']['candidate_count']} (τ, Δτ) pairs × {summary['grid']['rho_count']} ρ values per window.",
        f"- κ long 0.8 chip: {long_kappa['kappa']:.8g}; short-delay: {short_kappa['kappa']:.8g}.",
        f"- 0.8-chip gate clean active rate: {summary['gate_rates']['long_0p8/clean']['active_rate']:.3f}.",
        "- 0.8-chip counterfactual coverage: not scoreable from retained data (no same-second paired clean RINEX).",
        f"- Short-delay interval coverage: {summary['short_delay']['coverage_rate']:.3f}; median width {summary['short_delay']['median_width_m']:.3f} m.",
        f"- Stage-2 decision: **{summary['stage2_decision']}** under the frozen ≥90% rule.",
        "",
        "Every validation residual surface is retained as a compressed NPZ under `residual_surfaces/`; exact per-window endpoints are in `per_window_bias_intervals.csv`.",
    ]
    (RESULTS / "REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": summary["status"],
        "long_kappa": long_kappa["kappa"],
        "short_kappa": short_kappa["kappa"],
        "short_coverage": short_coverage,
        "short_median_width_m": summary["short_delay"]["median_width_m"],
        "stage2_decision": summary["stage2_decision"],
    }, indent=2), flush=True)


if __name__ == "__main__":
    main()
