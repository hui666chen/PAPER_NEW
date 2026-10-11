"""Cross-fitted, observation-only confidence for the two-path profile winner."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import ndtri
from scipy.stats import chi2, qmc, spearmanr

import compare_tau_pseudorange_proxy as tau_compare
import monte_carlo_resolution_baseline as mc
import run_confirmatory_outlier_mc as confirm_mc
import run_preregistered_outlier_mc as frozen_mc
from b_joint_blind_audit import SEPARATION, TAU
from blind_response_experiment import write_new
from run_blind_transfer import ROOT


DEFAULT_OUT = ROOT / "crossfit_basin_pilot_3x50"


class CrossfitBasinScorer:
    """Score observations without receiving simulation truth."""

    def __init__(self, dictionary, indices, inverse, competitors=8,
                 sobol_points=4096, sobol_scrambles=2):
        self.dictionary = np.asarray(dictionary, complex)
        self.indices = np.asarray(indices, int)
        self.inverse = np.asarray(inverse, complex)
        self.competitors = int(competitors)
        self.sobol_points = int(sobol_points)
        self.sobol_scrambles = int(sobol_scrambles)
        power = int(round(np.log2(self.sobol_points)))
        if 2 ** power != self.sobol_points:
            raise ValueError("sobol_points must be a power of two")
        self.sobol_power = power
        self.candidate_tau = np.tile(TAU, len(SEPARATION))
        self.candidate_separation = np.repeat(SEPARATION, len(TAU))

    def profile(self, correlations, weights):
        moment = np.einsum(
            "e,eq,es->qs", weights, correlations.conj(), correlations,
            optimize=True,
        )[None, :, :]
        costs = mc.candidate_costs(
            moment, np.zeros(1), self.indices, self.inverse
        )[0]
        return costs, int(np.argmin(costs))

    def discover(self, costs, anchor):
        order = np.argsort(costs, kind="stable")
        selected = []
        anchor_tau = self.candidate_tau[anchor]
        for candidate in order:
            candidate = int(candidate)
            if candidate == anchor:
                continue
            if abs(self.candidate_tau[candidate] - anchor_tau) <= .05 + 1e-12:
                continue
            if any(np.hypot(
                self.candidate_tau[candidate] - self.candidate_tau[old],
                self.candidate_separation[candidate]
                - self.candidate_separation[old],
            ) <= .05 + 1e-12 for old in selected):
                continue
            selected.append(candidate)
            if len(selected) == self.competitors:
                break
        return np.asarray(selected, int)

    def explained(self, correlations, candidate):
        cols = self.indices[candidate]
        values = correlations[:, cols]
        return np.einsum(
            "ei,ij,ej->e", values.conj(), self.inverse[candidate], values,
            optimize=True,
        ).real

    def heldout_probability(self, correlations, weights, anchor, competitors,
                            seed):
        if len(competitors) == 0 or len(correlations) < 3:
            return 1.0, [1.0], False
        anchor_explained = self.explained(correlations, anchor)
        differences = np.column_stack([
            weights * (
                anchor_explained - self.explained(correlations, int(candidate))
            )
            for candidate in competitors
        ])
        mean = differences.mean(axis=0)
        covariance = np.atleast_2d(np.cov(differences, rowvar=False, ddof=1))
        if not np.all(np.isfinite(mean)) or not np.all(np.isfinite(covariance)):
            return 1.0, [1.0], False
        covariance = (covariance + covariance.T) / 2 / len(differences)
        eigenvalues, eigenvectors = np.linalg.eigh(covariance)
        floor = max(float(eigenvalues[-1]) * 1e-10, 1e-18)
        factor = eigenvectors * np.sqrt(np.maximum(eigenvalues, floor))
        df = len(differences) - 1
        probabilities = []
        for scramble in range(self.sobol_scrambles):
            sampler = qmc.Sobol(
                d=len(mean) + 1, scramble=True, seed=int(seed + scramble)
            )
            uniform = np.clip(
                sampler.random_base2(self.sobol_power), 1e-12, 1 - 1e-12
            )
            normal = ndtri(uniform[:, :len(mean)]) @ factor.T
            scale = np.sqrt(df / chi2.ppf(uniform[:, -1], df))
            draw = mean[None, :] + normal * scale[:, None]
            probabilities.append(float(np.mean(np.min(draw, axis=1) < 0)))
        return float(np.mean(probabilities)), probabilities, True

    def direction(self, train_corr, train_weights, test_corr, test_weights,
                  seed):
        costs, anchor = self.profile(train_corr, train_weights)
        competitors = self.discover(costs, anchor)
        probability, scrambles, valid = self.heldout_probability(
            test_corr, test_weights, anchor, competitors, seed
        )
        return dict(
            probability=probability, valid=bool(valid), anchor=anchor,
            anchor_tau=float(self.candidate_tau[anchor]),
            anchor_separation=float(self.candidate_separation[anchor]),
            anchor_tau_boundary=bool(
                abs(self.candidate_tau[anchor] - TAU[0]) < 1e-12
                or abs(self.candidate_tau[anchor] - TAU[-1]) < 1e-12
            ),
            competitors=[int(v) for v in competitors],
            scramble_values=[float(v) for v in scrambles],
        )

    def score(self, epoch_correlations, weights, seed=1):
        correlations = np.asarray(epoch_correlations, complex)
        weights = np.asarray(weights, float)
        if len(correlations) != len(weights) or len(correlations) < 6:
            raise ValueError("correlations and weights must have equal length >= 6")
        index_a = np.arange(0, len(weights), 2)
        index_b = np.arange(1, len(weights), 2)
        full_costs, full_anchor = self.profile(correlations, weights)
        del full_costs
        a_to_b = self.direction(
            correlations[index_a], weights[index_a],
            correlations[index_b], weights[index_b], seed,
        )
        b_to_a = self.direction(
            correlations[index_b], weights[index_b],
            correlations[index_a], weights[index_a], seed + 100,
        )
        full_tau = float(self.candidate_tau[full_anchor])
        full_boundary = bool(
            abs(full_tau - TAU[0]) < 1e-12 or abs(full_tau - TAU[-1]) < 1e-12
        )
        consistency_failure = bool(
            full_boundary or a_to_b["anchor_tau_boundary"]
            or b_to_a["anchor_tau_boundary"]
            or not a_to_b["valid"] or not b_to_a["valid"]
            or abs(a_to_b["anchor_tau"] - full_tau) > .05 + 1e-12
            or abs(b_to_a["anchor_tau"] - full_tau) > .05 + 1e-12
            or abs(a_to_b["anchor_tau"] - b_to_a["anchor_tau"]) > .05 + 1e-12
        )
        combined = 1.0 if consistency_failure else max(
            a_to_b["probability"], b_to_a["probability"]
        )
        return dict(
            P_cross=float(combined), full_anchor=int(full_anchor),
            full_tau=full_tau, full_tau_boundary=full_boundary,
            consistency_failure=consistency_failure,
            P_A_to_B=float(a_to_b["probability"]),
            P_B_to_A=float(b_to_a["probability"]),
            tau_A=float(a_to_b["anchor_tau"]),
            tau_B=float(b_to_a["anchor_tau"]),
            A_to_B=a_to_b, B_to_A=b_to_a,
        )


def setup():
    (_, eigenvalues, eigenvectors, delay_grid, dictionary,
     indices, inverse, gain_pool) = frozen_mc.setup()
    return eigenvalues, eigenvectors, delay_grid, dictionary, indices, inverse, gain_pool


def run_cell(row, trials, scorer, model):
    eigenvalues, eigenvectors, delay_grid, dictionary, indices, inverse, gain_pool = model
    delay, amp, epochs = float(row["delay"]), float(row["amp"]), int(row["N_epoch"])
    trial_indices = np.arange(int(trials))
    noise_corr, noise_energy, weights = confirm_mc.synthesize_trials(
        trial_indices, delay, amp, epochs, eigenvalues, eigenvectors,
        dictionary, gain_pool,
    )
    base = mc.noise_moments(noise_corr, noise_energy, weights, epochs)
    direct_index = int(round((0 - delay_grid[0]) / .01))
    reflected_index = int(round((delay - delay_grid[0]) / .01))
    rho = amp * np.exp(1j * mc.PHASE_RAD)
    h = dictionary[:, direct_index] + rho * dictionary[:, reflected_index]
    hcorr = dictionary.conj().T @ h
    henergy = float(np.vdot(h, h).real)
    moment, _ = mc.moments(
        base, hcorr, henergy, [(1, direct_index), (rho, reflected_index)]
    )
    direct, single_indices, single_inverse, direct_norm2 = (
        tau_compare.build_single_path_models(delay_grid, dictionary)
    )
    single_tau, _ = tau_compare.baseline_tau(
        moment, direct, single_indices, single_inverse, direct_norm2
    )
    epoch_correlations = noise_corr[:, :epochs] + hcorr[None, None, :]
    records = []
    for trial in range(int(trials)):
        score = scorer.score(
            epoch_correlations[trial], weights[trial, :epochs],
            seed=20264000 + trial,
        )
        score.update(
            trial=int(trial), two_tau=float(score["full_tau"]),
            single_tau=float(single_tau[trial]),
            two_minus_single_squared_error=float(
                score["full_tau"] ** 2 - single_tau[trial] ** 2
            ),
            two_worse_than_single=bool(
                abs(score["full_tau"]) > abs(single_tau[trial]) + 1e-12
            ),
        )
        records.append(score)
    return dict(delay=delay, amp=amp, N_epoch=epochs, records=records)


def correlation(records):
    result = spearmanr(
        [r["P_cross"] for r in records],
        [r["two_minus_single_squared_error"] for r in records],
    )
    value = getattr(result, "statistic", result.correlation)
    return None if not np.isfinite(value) else float(value)


def summarize_cell(cell):
    records = cell["records"]
    false_trust = [
        r["two_worse_than_single"] and r["P_cross"] <= .20 + 1e-12
        for r in records
    ]
    return dict(
        delay=cell["delay"], amp=cell["amp"], N_epoch=cell["N_epoch"],
        trials=len(records), spearman=correlation(records),
        consistency_rate=float(np.mean([
            not r["consistency_failure"] for r in records
        ])),
        two_worse_rate=float(np.mean([
            r["two_worse_than_single"] for r in records
        ])),
        false_trust_rate=float(np.mean(false_trust)),
        low_risk_rate=float(np.mean([r["P_cross"] <= .20 for r in records])),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--max-cells", type=int, default=3)
    parser.add_argument("--trials", type=int, default=50)
    parser.add_argument("--competitors", type=int, default=8)
    parser.add_argument("--sobol-points", type=int, default=4096)
    parser.add_argument("--sobol-scrambles", type=int, default=2)
    args = parser.parse_args()
    args.output.mkdir(exist_ok=False)
    model = setup()
    scorer = CrossfitBasinScorer(
        model[3], model[4], model[5], competitors=args.competitors,
        sobol_points=args.sobol_points, sobol_scrambles=args.sobol_scrambles,
    )
    prediction_path = (
        ROOT / "outlier_theory_confirmatory_prediction_27_gaussian_frozen/predictions.json"
    )
    rows = json.loads(prediction_path.read_text())["results"][:args.max_cells]
    cells = []
    for number, row in enumerate(rows, start=1):
        cell = run_cell(row, args.trials, scorer, model)
        cells.append(cell)
        print(f"cell {number}/{len(rows)} {summarize_cell(cell)}", flush=True)
    flat = [record for cell in cells for record in cell["records"]]
    summaries = [summarize_cell(cell) for cell in cells]
    pooled_rho = correlation(flat)
    n20 = [row for row in summaries if row["N_epoch"] == 20]
    n20_supported = bool(n20) and all(
        row["consistency_rate"] >= .50 and row["false_trust_rate"] <= .10
        for row in n20
    )
    false_trust_rate = float(np.mean([
        r["two_worse_than_single"] and r["P_cross"] <= .20 + 1e-12
        for r in flat
    ]))
    checks = dict(
        pooled_spearman_ge_0p6405=bool(
            pooled_rho is not None and pooled_rho >= .6405
        ),
        every_cell_spearman_positive=bool(all(
            row["spearman"] is not None and row["spearman"] > 0
            for row in summaries
        )),
        N20_supported=n20_supported,
        pooled_false_trust_below_10_percent=bool(false_trust_rate < .10),
    )
    payload = dict(
        schema=1, method="cross-fitted observation-only basin confidence",
        truth_passed_to_score=False, endpoint_truth_used_after_scoring=True,
        split="zero-based even epochs A, odd epochs B",
        combination="max(P_A_to_B, P_B_to_A); force 1 on consistency failure",
        competitors=args.competitors,
        sobol=dict(points=args.sobol_points, scrambles=args.sobol_scrambles),
        pooled_spearman=pooled_rho, pooled_false_trust_rate=false_trust_rate,
        cell_summaries=summaries, expansion_checks=checks,
        expand_to_54=bool(all(checks.values())), cells=cells,
    )
    write_new(args.output / "results.json", payload)
    print(json.dumps({k: v for k, v in payload.items() if k != "cells"}, indent=2))
    print(f"COMPLETE={args.output}")


if __name__ == "__main__":
    main()
