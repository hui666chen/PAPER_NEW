#!/usr/bin/env python3
"""Chronological hold-out validation for the two multisatellite models."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

from compare_multisat_kalman import (
    CHIP_METERS,
    common_design,
    grouped,
    initialize_common,
    load_observations,
    load_truth,
)


def common_step(
    state: np.ndarray,
    covariance: np.ndarray,
    previous_window: int,
    window: int,
    group: list[dict[str, float | int | str]],
    truth: dict[int, np.ndarray],
    q_pos_m: float,
    q_clock_m: float,
) -> tuple[float, np.ndarray, np.ndarray]:
    steps = max(int(window) - int(previous_window), 1)
    process = np.diag([q_pos_m * q_pos_m * steps] * 3 + [q_clock_m * q_clock_m * steps])
    predicted_covariance = covariance + process
    h, y, r = common_design(group, truth)
    innovation = y - h @ state
    innovation_covariance = h @ predicted_covariance @ h.T + np.diag(r)
    sign, logdet = np.linalg.slogdet(innovation_covariance)
    if sign <= 0.0 or not math.isfinite(logdet):
        return float("inf"), state, covariance
    try:
        solved = np.linalg.solve(innovation_covariance, innovation)
        gain = np.linalg.solve(innovation_covariance, h @ predicted_covariance).T
    except np.linalg.LinAlgError:
        return float("inf"), state, covariance
    n = len(y)
    nll = 0.5 * (n * math.log(2.0 * math.pi) + logdet + float(innovation @ solved))
    updated_state = state + gain @ innovation
    identity = np.eye(4)
    residual_map = identity - gain @ h
    updated_covariance = residual_map @ predicted_covariance @ residual_map.T + gain @ np.diag(r) @ gain.T
    updated_covariance = 0.5 * (updated_covariance + updated_covariance.T)
    return float(nll), updated_state, updated_covariance


def independent_step(
    state: np.ndarray,
    covariance: np.ndarray,
    previous_window: int,
    window: int,
    group: list[dict[str, float | int | str]],
    prns: list[int],
    q_m: float,
    rho: float,
) -> tuple[float, np.ndarray, np.ndarray]:
    index = {prn: i for i, prn in enumerate(prns)}
    steps = max(int(window) - int(previous_window), 1)
    transition = rho**steps
    process = q_m * q_m * (1.0 - transition * transition) if abs(rho) < 1.0 else q_m * q_m * steps
    predicted_state = transition * state
    predicted_covariance = transition * transition * covariance + process
    nll = 0.0
    for row in group:
        i = index[int(row["prn"])]
        measurement = float(row["y"]) * CHIP_METERS
        measurement_variance = max(float(row["sigma"]) * CHIP_METERS, 1.0e-3) ** 2
        innovation = measurement - predicted_state[i]
        innovation_variance = predicted_covariance[i] + measurement_variance
        if innovation_variance <= 0.0 or not math.isfinite(innovation_variance):
            return float("inf"), state, covariance
        nll += 0.5 * (
            math.log(2.0 * math.pi * innovation_variance)
            + innovation * innovation / innovation_variance
        )
        gain = predicted_covariance[i] / innovation_variance
        predicted_state[i] += gain * innovation
        predicted_covariance[i] = (1.0 - gain) * predicted_covariance[i]
    return float(nll), predicted_state, predicted_covariance


def common_train_test(
    train: list[tuple[int, list[dict[str, float | int | str]]]],
    test: list[tuple[int, list[dict[str, float | int | str]]]],
    truth: dict[int, np.ndarray],
    q_pos_m: float,
    q_clock_m: float,
    prior_std_m: float,
) -> tuple[float, float]:
    state, covariance = initialize_common(train[0][1], truth, prior_std_m)
    previous_window = train[0][0]
    train_nll = 0.0
    for window, group in train[1:]:
        nll, state, covariance = common_step(
            state, covariance, previous_window, window, group, truth, q_pos_m, q_clock_m
        )
        train_nll += nll
        previous_window = window
    test_nll = 0.0
    for window, group in test:
        nll, state, covariance = common_step(
            state, covariance, previous_window, window, group, truth, q_pos_m, q_clock_m
        )
        test_nll += nll
        previous_window = window
    return float(train_nll), float(test_nll)


def independent_train_test(
    train: list[tuple[int, list[dict[str, float | int | str]]]],
    test: list[tuple[int, list[dict[str, float | int | str]]]],
    prns: list[int],
    q_m: float,
    rho: float,
    prior_std_m: float,
) -> tuple[float, float]:
    index = {prn: i for i, prn in enumerate(prns)}
    state = np.zeros(len(prns), dtype=np.float64)
    covariance = np.full(len(prns), prior_std_m**2, dtype=np.float64)
    for row in train[0][1]:
        i = index[int(row["prn"])]
        state[i] = float(row["y"]) * CHIP_METERS
        covariance[i] = max(float(row["sigma"]) * CHIP_METERS, 1.0e-3) ** 2
    previous_window = train[0][0]
    train_nll = 0.0
    for window, group in train[1:]:
        nll, state, covariance = independent_step(
            state, covariance, previous_window, window, group, prns, q_m, rho
        )
        train_nll += nll
        previous_window = window
    test_nll = 0.0
    for window, group in test:
        nll, state, covariance = independent_step(
            state, covariance, previous_window, window, group, prns, q_m, rho
        )
        test_nll += nll
        previous_window = window
    return float(train_nll), float(test_nll)


def search(
    train: list[tuple[int, list[dict[str, float | int | str]]]],
    test: list[tuple[int, list[dict[str, float | int | str]]]],
    truth: dict[int, np.ndarray],
    prns: list[int],
    args: argparse.Namespace,
) -> dict[str, object]:
    q_common = np.geomspace(args.q_min_m, args.q_max_m, args.grid_size)
    best_common = None
    for q_pos in q_common:
        for q_clock in q_common:
            train_nll, test_nll = common_train_test(
                train, test, truth, float(q_pos), float(q_clock), args.prior_std_m
            )
            if best_common is None or train_nll < best_common[0]:
                best_common = (train_nll, test_nll, float(q_pos), float(q_clock))
    assert best_common is not None

    q_gm = np.geomspace(args.gm_q_min_m, args.gm_q_max_m, args.grid_size)
    rhos = np.linspace(args.rho_min, args.rho_max, args.rho_grid_size)
    best_gm = None
    for q in q_gm:
        for rho in rhos:
            train_nll, test_nll = independent_train_test(
                train, test, prns, float(q), float(rho), args.prior_std_m
            )
            if best_gm is None or train_nll < best_gm[0]:
                best_gm = (train_nll, test_nll, float(q), float(rho))
    assert best_gm is not None

    common_test_ll = -best_common[1]
    gm_test_ll = -best_gm[1]
    return {
        "train_windows": [window for window, _ in train],
        "test_windows": [window for window, _ in test],
        "n_train_observations": sum(len(group) for _, group in train),
        "n_test_observations": sum(len(group) for _, group in test),
        "common_state": {
            "train_log_predictive_likelihood": -best_common[0],
            "test_log_predictive_likelihood": common_test_ll,
            "q_position_m": best_common[2],
            "q_clock_m": best_common[3],
        },
        "independent_gauss_markov": {
            "train_log_predictive_likelihood": -best_gm[0],
            "test_log_predictive_likelihood": gm_test_ll,
            "q_m": best_gm[2],
            "rho": best_gm[3],
        },
        "test_loglik_difference_common_minus_gm": common_test_ll - gm_test_ll,
        "decision": "COMMON_STATE" if common_test_ll > gm_test_ll else "INDEPENDENT_GAUSS_MARKOV",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--truth", type=Path, required=True)
    parser.add_argument("--estimates", type=Path, required=True)
    parser.add_argument("--scenario", choices=("spoof", "multipath"), required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--prior-std-m", type=float, default=2000.0)
    parser.add_argument("--q-min-m", type=float, default=1.0e-5)
    parser.add_argument("--q-max-m", type=float, default=1000.0)
    parser.add_argument("--gm-q-min-m", type=float, default=1.0e-4)
    parser.add_argument("--gm-q-max-m", type=float, default=1000.0)
    parser.add_argument("--rho-min", type=float, default=0.50)
    parser.add_argument("--rho-max", type=float, default=0.9999)
    parser.add_argument("--rho-grid-size", type=int, default=25)
    parser.add_argument("--grid-size", type=int, default=13)
    parser.add_argument("--split-index", type=int, default=None)
    args = parser.parse_args()

    truth = load_truth(args.truth)
    rows = load_observations(args.estimates, truth, args.scenario)
    groups = grouped(rows)
    if len(groups) < 4:
        raise ValueError("at least four windows are required for hold-out validation")
    split_index = args.split_index if args.split_index is not None else len(groups) // 2
    split_index = max(2, min(split_index, len(groups) - 2))
    train = groups[:split_index]
    test = groups[split_index:]
    prns = sorted({int(row["prn"]) for row in rows})
    result = search(train, test, truth, prns, args)
    result.update(
        {
            "scenario": args.scenario,
            "prns": prns,
            "split_index": split_index,
            "n_windows": len(groups),
            "chip_meters": CHIP_METERS,
            "parameter_search": {
                "q_position_m": [args.q_min_m, args.q_max_m],
                "q_gm_m": [args.gm_q_min_m, args.gm_q_max_m],
                "rho": [args.rho_min, args.rho_max],
            },
        }
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(f"scenario={args.scenario}")
    print(f"train_windows={result['train_windows']} test_windows={result['test_windows']}")
    print(f"n_train_observations={result['n_train_observations']} n_test_observations={result['n_test_observations']}")
    print(f"common_test_log_predictive_likelihood={result['common_state']['test_log_predictive_likelihood']:.6f}")
    print(f"independent_gm_test_log_predictive_likelihood={result['independent_gauss_markov']['test_log_predictive_likelihood']:.6f}")
    print(f"test_loglik_difference_common_minus_gm={result['test_loglik_difference_common_minus_gm']:.6f}")
    print(f"decision={result['decision']}")
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
