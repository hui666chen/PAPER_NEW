# Outlier-theory follow-up report

Date: 2026-09-25

## Frozen 27-cell comparison

The same Monte Carlo outcomes and the same preregistered event definitions were used for both predictors. No B-batch sealed truth was accessed.

| Predictor | Local Wilson pass | Global Wilson pass | Local Spearman (95% bootstrap CI) | Global Spearman (95% bootstrap CI) |
|---|---:|---:|---:|---:|
| Gaussian moment approximation | 25/27 | 26/27 | 0.993 (0.963, 1.000) | 0.986 (0.952, 0.997) |
| Saddlepoint-calibrated approximation | 24/27 | 25/27 | 0.994 (0.966, 1.000) | 0.984 (0.942, 0.996) |

The saddlepoint implementation uses a Lugannani--Rice CDF for every indefinite noncentral quadratic-form marginal and retains the previously frozen Gaussian copula for joint events. It is therefore a marginal saddlepoint correction, not a full multivariate saddlepoint approximation.

Across all 27 cells, the low-rank quadratic representation reproduced the existing analytic means and variances with maximum relative discrepancies of 1.81e-8 and 1.89e-7, respectively. No marginal required a Gaussian fallback. The largest absolute change from the moment predictor was 0.00870 for the local probability and 0.01045 for the global probability.

Conclusion: marginal non-Gaussianity is not the dominant source of prediction error on the frozen boundary map. The saddlepoint correction does not remove the preregistered trigger and should be reported as a negative ablation rather than selected as the main model.

This failure to remove the trigger is logically consistent with the original diagnosis. The active stratified-Spearman trigger came from only three cells in the `<1e-4` layer, whose perfect rank ordering produced a degenerate bootstrap interval. A tail-probability approximation cannot repair that sample-size defect. Its failure to do so therefore supports, rather than contradicts, the diagnosis that the trigger was caused by sparse stratification.

### Changed-cell implementation audit

The two cells that changed from Wilson pass to fail were audited separately. For the four most competitive candidates in each cell, Lugannani--Rice marginal CDFs were compared with direct Gil--Pelaez inversion of the exact empirical-weight-mixture characteristic function. Across the eight candidates, the Lugannani--Rice absolute error was `2.2e-6` to `2.5e-4`, and every result was more accurate than its Gaussian marginal. Direct inversion was unchanged between integration limits 160 and 320. Thus saddlepoint root finding was stable.

The degraded joint-event score must instead be attributed to the remaining joint approximation: corrected marginal CDFs were coupled with the old Gaussian copula. This audit validates the marginal saddlepoint implementation but does not test a full multivariate saddlepoint approximation. Keeping the Gaussian moment predictor is justified as the better calibrated complete predictor on the frozen 27 cells; it is not a claim that a full multivariate saddlepoint method has been disproved.

### Empirical pair-copula audit

The exact frozen 200-trial base stream and 800-trial supplement stream were replayed for both changed-status cells. For each cell, all six pairs among the four candidates chosen analytically before this replay were examined. At lower-tail quantiles `u = 0.02, 0.05, 0.10, 0.20`, all 48 matched-correlation Gaussian-copula probabilities lay inside the empirical Wilson 95% intervals. The same was true for all 12 candidate-pair zero-threshold probabilities. The maximum absolute finite-quantile discrepancy was 0.00988 and the median was 0.00140.

Consequently, the available data do not support the stronger claim that pairwise Gaussian-copula tail dependence is structurally wrong. The remaining discrepancy may be a high-dimensional extreme-value/competition effect accumulated over hundreds of local and remote candidates. Although the exact joint CGF of a fixed set of correlated Hermitian quadratic forms is analytic (and remains available under the empirical epoch-weight mixture), evaluating the complete event is still a high-dimensional orthant/minimum problem. A reduced active-set multivariate inversion or saddlepoint calculation is technically natural, but it should not be launched solely on a pairwise-tail-mismatch claim that this audit did not confirm.

## Stratified supplementary design

An initial 18-condition design was produced without using Monte Carlo outcomes. Before any of those cells were simulated or scored, it was superseded by a 27-condition confirmatory design with nine cells in each normalized remote-cost-gap tier (`<1e-4`, `1e-4~3e-4`, and `>3e-4`). The design deliberately fills the previously empty middle tier and balances epoch counts within every tier. The exact cells and deterministic selection audit are stored in `gps_l1_prn14_blind_transfer/outlier_theory_stratified_design_v2/design.json` and frozen in `PREREGISTER_outlier_theory_confirmatory.md`.

The middle-tier candidate pool is smaller (26 versus 132 and 164), but the selected middle-tier cost gaps are not more widely spaced. Their median adjacent gap is `2.05e-5`, compared with `1.86e-4` and `2.03e-4` in the low and high tiers; the median adjacent distance in `(delay, amplitude)` is 0.137, compared with 0.164 and 0.214. The middle tier is therefore thinner in total range, not sparser after selection. This distinction should be stated when comparing trends across tiers.

The frozen 27-cell rescore still has zero middle-tier cells. Therefore the stratified trigger cannot become fully evaluable until these supplementary cells are simulated and scored; changing the approximation alone cannot repair that design gap.

## Remote-anchor full-enumeration audit

For `(delay=0.60, amp=0.10, N=20)`, 200 matched trials gave global-event rates 0.700 with all 8,220 remote candidates and 0.695 with the 505-candidate anchor approximation. Winner agreement was 99.5%, with one anchor false negative and no false positives.

For `(delay=0.80, amp=0.20, N=100)`, both methods gave zero global events in 200 matched trials and 100% agreement.

Conclusion: the row-anchor reduction is not responsible for the observed 27-cell calibration behavior at the two audited operating points.

## Delay-controlled rank analysis

After rank-transforming the variables and regressing both prediction and Monte Carlo frequency on injected delay rank, the partial Spearman correlations were:

- Local: 0.9939, bootstrap 95% CI (0.9692, 0.9996)
- Global: 0.9860, bootstrap 95% CI (0.9509, 0.9966)

Thus the strong predictor--outcome ordering is not merely induced by injected delay. This audit controls delay only; amplitude ratio and epoch count remain part of the designed operating grid.

## Decision

Keep the simpler Gaussian moment predictor as the current primary semi-analytic model. Run the new 27-cell independent confirmatory experiment next, then repeat the frozen local/global Wilson and tiered-Spearman scoring and add the preregistered exact permutation audit. The central unresolved issue is coverage of the cost-gap strata, especially the middle tier, rather than the lack of a saddlepoint marginal correction.

## Independent 27-cell confirmation

The Gaussian-moment predictions were frozen before outcome generation. Monte Carlo seeds were independently derived for every `(cell, trial)` pair from the frozen prediction hash. Twenty-four cells received 200 trials and the three prospectively selected cells with `P_global_pred < 0.02` received 1,000 trials directly; no interim outcome was used to choose the trial count.

The confirmatory Wilson pass counts were 23/27 (85.2%) for local events and 26/27 (96.3%) for global events. Overall Spearman correlations were 0.960 for local probabilities and 0.983 for global probabilities. The three frozen upgrade triggers were all false.

Every tier contained nine cells and all `9! = 362,880` label permutations were enumerated for the exact Spearman audit. Tier-level local Spearman values were 0.933, 0.667, and 0.996 from low to high cost gap, with exact two-sided p-values 0.00075, 0.0589, and 0.000011. Global values were 0.933, 0.833, and 0.983, with p-values 0.00075, 0.00827, and 0.000050. The middle-tier local association is the weakest and its exact p-value is slightly above 0.05, but its bootstrap interval contains the overall coefficient and its point estimate exceeds the frozen half-overall threshold; it therefore does not trigger an upgrade under the preregistered rule.

This confirmation resolves the sparse-tier artifact that triggered the earlier 27-cell audit. The data support retaining the Gaussian moment predictor as the primary semi-analytic outlier-probability model, with the saddlepoint and copula analyses reported as diagnostic ablations.

## Direct-path delay and pseudorange-proxy endpoint

The MC winner coordinate `tau` is the estimated direct-path delay, with injected truth `tau_0 = 0`. GPS L1 C/A conversion uses `c / 1.023 MHz = 293.0523 m/chip`; the 0.01-chip estimator grid is 2.9305 m. The sealed B truth contains relative path delay, amplitude, and phase but no absolute satellite--receiver geometric range, so the reported metre endpoint is direct-path delay error converted to metres, not a complete absolute pseudorange or PVT error.

Two baselines were replayed on the exact same random trials: a single-path global NLS with the same response dictionary and nuisance columns, and a normalized noncoherent correlation-energy peak proxy. Across the old 27 cells, the two-path estimator beat single-path NLS in 19 cells and reduced mean per-cell tau RMSE from 24.19 m to 15.82 m. Across the deliberately difficult new 27 cells, it beat single-path NLS in only 7 cells and increased mean per-cell RMSE from 13.19 m to 35.25 m. Thus unconditional two-path correction is not deployable.

The 35.25 m figure is the arithmetic mean of the 27 cell-wise RMSE values,
`mean_j sqrt(mean_r(e_jr^2))`. It is not the same summary as the 41.53 m
equal-cell aggregate RMSE used below, `sqrt(mean_j mean_r(e_jr^2))`. Both give
every cell equal weight, but the latter squares before averaging and therefore
places more influence on high-error cells. The equal-cell aggregate RMSE is the
primary cross-batch endpoint; mean per-cell RMSE is retained only as a
descriptive cell-level summary.

The frozen global-outlier probability is strongly associated with the relative value of the two-path estimate. Spearman correlation between predicted risk and `(two-path RMSE - single-path RMSE)` is 0.957 in the old 27 cells and 0.915 in the new 27 cells. On a prelisted threshold grid, selecting the threshold on the old 27 cells alone gives `P_global <= 0.20`; applying this unchanged threshold to the new 27 cells reduces equal-cell aggregate tau RMSE from 14.15 m (always single path) and 41.53 m (always two path) to 11.40 m. A cell-only paired bootstrap gives a gate-minus-single 95% interval of `[-4.89, -0.89] m`. A two-level paired bootstrap, which first resamples cells and then resamples matched trials within every selected cell, gives `[-4.95, -0.84] m`; its gate-RMSE interval is `[9.57, 13.16] m`. The added within-cell Monte Carlo uncertainty therefore does not change this retrospective comparison.

This threshold analysis is a post-hoc old-to-new train/test audit, not a preregistered primary endpoint. More importantly, the present `P_global` is an oracle simulation quantity: its signal mean, reference candidate, and local/remote candidate sets are constructed from the injected true delay, amplitude ratio, and true `(tau_0, delta_tau_0)`. It therefore cannot yet be evaluated by a receiver without truth. The result establishes that this theoretical risk label is informative; it does not establish a deployable gate. Deployment requires an observation-only replacement, such as the fitted-profile competing-basin gap or a plug-in parametric-bootstrap basin-switch probability, calibrated on independent data and frozen before the subsequent receiver/IQ pseudorange or PVT evaluation.

An observation-only winner-centred analytic plug-in replacement was subsequently implemented and tested on 150 trials from three difficult confirmatory cells. It reproduced the saved estimator winners exactly, but its trial-level risk/loss Spearman correlation was -0.301; within the three cells the correlations ranged from -0.674 to -0.705. Active-set and Sobol-resolution sensitivity checks did not change the direction. The failure is structural: when the observed winner is an erroneous basin, fitting the data-generating mean at that same basin can assign high conditional stability to the wrong solution. Consequently, neither this naive plug-in score nor a parametric bootstrap centred on the same winner is accepted as the deployable gate. The next candidate must break selection reuse, for example through cross-fitted basin discovery/evaluation or an epoch/block bootstrap that reselects the winner.

A deterministic two-fold cross-fitted basin score was then frozen and tested
on the same three-by-50 stress pilot. It changed all three within-cell
risk/loss correlations from negative to positive (0.504, 0.728, and 0.449) and
produced no false-trust trial at the prospectively specified risk threshold
0.20. However, pooled Spearman was 0.558, below the frozen expansion criterion
0.6405, and only 4--12% of trials were assigned risk at or below 0.20 because
split/full winner consistency was only 40--52%. The N=20 cell narrowly passed
its separate support rule with 52% consistency and zero false trust. Under the
all-conditions-required protocol, the method is not expanded to 54 cells and
is characterized as safe but insufficiently informative. No split,
combination, threshold, or consistency rule was retuned after this result.
