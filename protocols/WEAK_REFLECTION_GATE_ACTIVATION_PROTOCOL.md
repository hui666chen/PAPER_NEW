# Frozen weak-reflection gate-activation protocol

Status: exploratory synthetic-observation study designed before evaluating the
amplitude-ratio 0.01 outcomes. It uses no new IQ capture.

## Why simulation, not the oracle formula

The existing analytic `P_global` is truth-centred and cannot be interpreted as
the deployment activation rate of the observation-only cross-fitted gate.
Therefore activation is evaluated by passing synthetic noisy observations
through the unchanged `P_cross` implementation. Oracle `P_global` is not used
to make the gate decision.

## Frozen grid

- amplitude ratio: exactly 0.01, corresponding to -40 dB power;
- reflection phase: the existing frozen phase in
  `monte_carlo_resolution_baseline.py`;
- delays: 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.60, 0.80 chip;
- epoch counts: 20, 50, 100;
- trials: first 50 deterministic confirmatory-seed trials per cell;
- total: 24 cells and 1,200 trials;
- response dictionary, covariance, epoch-weight distribution, profile search,
  split, active set, Sobol integration, consistency fallback, and gate
  threshold `P_cross <= 0.20`: unchanged.

Truth parameters are used only by the observation generator and endpoint
scoring. They are not passed to `CrossfitBasinScorer.score`.

## Primary endpoint and frozen labels

Primary endpoint: proportion of trials in which the gate selects the two-path
estimate.

- **Nearly always fallback:** pooled activation <= 5% and every individual
  cell activation <= 10%.
- **Mostly fallback:** pooled activation <= 10% but the previous stricter rule
  fails.
- **Non-negligible activation:** pooled activation > 10%.

Also report activation by delay and epoch count, split/full consistency,
selected-better/tie/worse counts, and tau-proxy RMSE.

## Safety criterion

A false-trust trial is one in which the gate activates and the two-path
absolute tau error is strictly worse than the single-path error. The weak-case
safety-net criterion passes if the pooled false-trust rate is <= 1% of all
1,200 trials. This safety result is descriptive and does not supersede the
earlier failed 54-cell expansion criterion.

## Implication for the future simulator matrix

If the result is nearly always or mostly fallback, the company multi-satellite
matrix should contain both:

1. amplitude 0.01 as the realistic urban-canyon safety case;
2. separately labelled strengthened amplitudes 0.05 and 0.10 as correction-
   benefit cases.

The strengthened cases must not be presented as the original region-C model.
