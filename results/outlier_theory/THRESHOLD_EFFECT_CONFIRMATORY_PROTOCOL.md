# Frozen confirmatory protocol: correlator-domain threshold-effect model

Frozen on 2026-09-25 before generating the confirmatory Monte Carlo outcomes.
This protocol concerns the A-derived synthetic experiment only. B truth remains
sealed and is not an input.

## Development and confirmation split

The existing 96-cell/100-repeat baseline, the selected 27-cell/200-repeat
extension, and the six-cell semi-analytic audit are development evidence. They
must not be presented as a pre-registered confirmation of the final predictor.

Confirmation uses 200 newly generated repeats in each of the same 27 design
cells, with RNG seed `2026092501`. The response, covariance, estimator grid,
phase, epoch-weight distribution, and cell list remain frozen. Common random
numbers are used across cells for paired comparisons; repeats within a cell
are independent. Cell-bootstrap intervals therefore describe variation across
the 27 fixed design points, not independent Monte Carlo noise across cells.

## Operational outcome classes

Let `e = abs(separation_hat - separation_true)` and use `1e-12` chip only to
protect exact grid thresholds from binary floating-point error.

- correct: `e <= 0.05 + 1e-12` chip;
- local-frontier outlier: `0.05 + 1e-12 < e <= 0.10 + 1e-12` chip;
- remote/global-operational outlier: `e > 0.10 + 1e-12` chip;
- total outlier: local-frontier plus remote/global-operational.

“Global-operational” is a scoring definition, not a claim that every estimate
beyond 0.10 chip belongs to a distinct topological basin. A separate watershed
analysis would be required for that stronger claim.

## Frozen predictor

For every candidate relative to the exact truth-grid candidate,

`D_i = C_i - C_0 = y^H (P_0 - P_i) y`.

The predictor uses the exact first two moments and cross-covariances of these
noncentral indefinite Hermitian quadratic forms under the frozen proper-complex
Gaussian covariance model. It does not fit coefficients to Monte Carlo rates.

Candidates are partitioned into correct, local-frontier, and remote regions.
The retained representation contains at most 30 correct, 60 local, and 60
remote candidates. Half of each budget is filled by the smallest standardized
mean gap; the remainder preserves deterministic grid-local minima and then the
next smallest gaps. Correlated Gaussian cost differences are integrated using
four independently scrambled Sobol runs of 32768 points each. The three class
probabilities are determined by which region has the smallest sampled cost.

Representation sensitivity doubles the local and remote budgets. An absolute
probability change above 0.01 in any class is a candidate-representation
failure and must be repaired before judging the quadratic-form tail model.

## Primary acceptance criteria

Local-frontier and remote/global-operational predictions are judged separately;
good total prediction cannot compensate for failure of either component.

For each of local, remote, and total probability:

1. at least 80% of the 27 predictions must lie inside the cellwise 95% Wilson
   interval from the 200 confirmatory repeats; and
2. mean absolute probability error must not exceed 0.04.

Both criteria must pass separately for local and remote. Total probability is
reported as a third primary check, not as a substitute.

## Rank correlation and approximation-stress diagnostics

Spearman correlation is computed separately for local, remote, and total rates,
with a 10,000-resample percentile interval obtained by resampling design cells.
It is used as a gate only when both observed and predicted rates contain at
least five distinct values. An eligible class triggers tail-model escalation
when its point estimate is below 0.60 or its interval lower endpoint is below
0.10.

The analytic stress score is the smallest pairwise standardized cost gap from
truth to any outlier candidate. Cells with score <= 1.645 form the threshold
subset; the others form the regular subset. If both subsets contain at least
six cells, a threshold-subset Spearman correlation at least 0.20 below the
regular-subset correlation, together with threshold-subset MAE above 0.04,
triggers tail-model escalation.

## Escalation rule

- Candidate-budget sensitivity failure: improve continuous-surface/basin
  representation first; do not blame the quadratic-form CDF.
- Otherwise, failure of either separate class criterion or an eligible rank/
  stress gate: replace the Gaussian moment tail with a noncentral indefinite
  quadratic-form saddlepoint or characteristic-function inversion.
- If all gates pass: retain the cheaper correlated Gaussian/Sobol predictor and
  validate it later against receiver-state-varying and hardware/simulator data.

## Frozen source hashes

- Noise covariance: `ba63b49d1829153eb4168a928acb9565a52c080a1b29141879bb4106f5564f73`
- A response waveform: `2884c11c7ef752a9f51a079a72def72bb96c78259a4863235281c31559568efe`
- Six-cell predictor at freeze: `cd71848b229f952cef70c7a05619fdc216a5a9bda149c5bb454c3d53da34ee04`
- Frozen-estimator MC implementation: `8e7e64870d9da71d4e2cf7ded5d4ef0de1da91ae16be4c980ff330247678d5fe`

