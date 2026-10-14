# Numerical implementation audit

Date: 2026-10-14 (Asia/Shanghai)

## Conclusion

The raw-Gram pseudoinverse does silently change the historical estimator's
refinement model: with `rcond=1e-12`, it drops the constant and linear nuisance
directions. But after replacing that hybrid refinement with a stable,
column-scaled four-column fit, the retained PRN30 batch kept the same coarse
basins, confidence probabilities, activation decisions, and persistence
decisions. The PRN30 0.10-amplitude PVT improvement also remained in a
same-observation re-run. Thus this audit does **not** support the claim that the
PRN30 positive result was merely a pseudoinverse artifact.

This conclusion is limited to one retained PRN30 realization. The original
PRN9 validation's full per-epoch response records were deleted after scoring,
so its ten-run result could not be reprocessed with the stable estimator and
remains un-audited by this exact test.

## What the old PVT code computes

The historical PRN9 and PRN30 pipelines use a four-column independent-path
design `[direct, reflected, constant, linear]` in the coarse candidate score.
That stage uses an unscaled `np.linalg.inv(Gram)`. Continuous tau refinement
then uses an unscaled `np.linalg.pinv(Gram, rcond=1e-12)`.

At representative PRN9/PRN30 candidates, the raw Gram's smallest-to-largest
eigenvalue ratio is about `3.1e-14` to `3.6e-14`. The refinement cutoff thus
keeps effective rank 2 of 4 and discards the constant and linear nuisance
directions. After column scaling, the four-column design condition number is
about 2.8--2.9 at 0.8 chip and about 13.9 at 0.1 chip. The observed truncation
therefore comes from scale disparity in the raw Gram, not a claim that the
column-normalized long-delay model is nearly singular.

The old estimator is hybrid: coarse ranking includes affine columns, while
continuous refinement effectively fits only the two path columns. This is a
real model change; it is not by itself proof that the final PVT result is
wrong.

## Shared-ratio synthetic check (separate estimator)

The three-column shared-ratio branch was tested separately on PRN9/PRN30 banks
at 0.8- and 0.1-chip separations. Its 200 trials per PRN/separation use
five-epoch block resampling of paired held-out clean-window affine
coefficients and residual vectors. These are empirical receiver/model
perturbations, not independent thermal-noise draws.

For this branch, stable B projects only onto the real-valued constant/linear
nuisance basis. The initial row-vector conjugation typo is algebraically
equivalent for that real basis, so the original shared-ratio numerical results
were unaffected; a corrected copy is retained with a projection assertion.
In the tested cases, old A and stable no-affine C agreed, while stable affine
B improved some separation estimates. This is a separate diagnostic and does
not establish performance of the four-column PVT estimator.

## Four-column known-truth check

Noiseless local-grid tests used the actual PRN9/PRN30 response banks, true
`tau=0`, separations `0.8` and `0.1` chip, `rho=0.34∠22.5°`, and either no
affine injection or constant/linear terms each set to 5% of path RMS. The
correct complex row projection was checked against `lstsq` before scoring.

- With no affine injection, all solvers recovered the on-grid truth in these
  local tests.
- With affine injection, stable four-column B recovered the truth in all
  tested cases. The historical pseudoinverse refinement A shifted the
  0.8-chip PRN9 case to `0.01/0.81` chip and shifted the 0.1-chip cases in
  several PRN/bank combinations; the stable two-column C showed the same
  shifted solutions as A.
- Short-delay profiles can have several nearly equivalent local candidates;
  one grid-step difference alone is not evidence of numerical failure.

This is a deterministic solver diagnostic, not a performance or gate test.

## Retained PRN30 data: 64-window re-score

All 64 retained windows were processed: 16 each for clean and amplitudes 0.01,
0.05, and 0.10. The response banks, input responses, delay/separation grids,
4096-point two-scramble cross-fit, and `P_cross <= 0.20` threshold were fixed.
Historical A reproduced its saved estimates exactly in all 64 windows.

| Scene | Cross-active windows A/B | Cross-active C | Max A-to-B `|Δtau_hat|` | Persistent 3-of-4 A/B |
|---|---:|---:|---:|---:|
| clean | 3 / 3 | 12 | 0.00526 chip | 0 / 0 |
| amplitude 0.01 | 4 / 4 | 12 | 0.00479 chip | 0 / 0 |
| amplitude 0.05 | 4 / 4 | 14 | 0.00664 chip | 0 / 0 |
| amplitude 0.10 | 7 / 7 | 12 | 0.00712 chip | 2 / 2 |

A/B `P_cross` values, coarse separations, and cross-gate decisions matched in
all 64 windows. The largest refined tau difference is `0.00712 chip`, about
`2.13 m` for an individual correction. C omits the affine nuisance model and
substantially changes the fitted basins and activation counts; C is therefore
a model-choice comparison, not just a better-conditioned implementation.
For each scene, the stable scaled-Gram coarse profile was also compared with a
direct complex-QR projection at three representative grid candidates; the
largest relative score discrepancy was below `9.5e-14`.

## PRN30 common-observation PVT sensitivity

The corrected B schedules were applied to copies of the retained common
single-path RINEX. Only G30 C1C changed; all other fields were byte-identical.
RTKLIB used the same executable, configuration, and navigation file. The
single baseline and both new solutions had identical 114 timestamp sets; each
scored scene contained all 16 expected epochs.

3D RMSE in metres (one retained realization):

| Scene | Single baseline | Historical A cross-only | Stable B cross-only | Historical A persistent | Stable B persistent |
|---|---:|---:|---:|---:|---:|
| clean | 2.1866 | 2.3800 | 2.2579 | 2.1866 | 2.1866 |
| amplitude 0.01 | 2.3827 | 2.3792 | 2.3578 | 2.3827 | 2.3827 |
| amplitude 0.05 | 2.8137 | 2.7106 | 2.7287 | 2.8137 | 2.8137 |
| amplitude 0.10 | 2.8540 | 2.3713 | 2.3480 | 2.7642 | 2.7656 |

The PRN30 amplitude-0.10 cross-only improvement remains (`2.8540` to
`2.3480 m`); the persistent result also remains essentially unchanged
(`2.7642` vs `2.7656 m`). In the clean scene, cross-only still worsens PVT,
while the existing 3-of-4 persistence rule applies no correction and equals
the single baseline. These are not independent-seed confidence intervals.

## Scope and remaining work

- The approximately 1.9 m clean-zero-point scatter from the separate
  stable-projection analysis is unaffected by the `pinv` issue and remains a
  valid descriptive number; this audit does not explain its cause.
- PRN30's retained-batch positive PVT result survives this numerical
  sensitivity check, but the result is one realization and does not establish
  seed-level robustness by itself.
- PRN9's full validation response matrices are unavailable; its reported
  ten-run result needs a new reproduction to be checked with stable
  four-column fitting.
- Whether the affine nuisance terms represent a physical offset, structured
  background, or response mismatch remains a model question. Their inclusion
  here is a controlled comparison, not a universal physical assertion.
- The short-delay covariance/FIM reconciliation is a separate follow-up and
  has not been completed here.

## Reproducible artifacts

- `synthetic_mc_200_corrected.json`: checked shared-ratio synthetic tests.
- `historical_four_column_noiseless_corrected.json`: four-column known-truth
  diagnostic.
- `prn30_retained_four_column_rescore_corrected.csv` and
  `prn30_retained_four_column_rescore_summary_corrected.json`: per-window
  PRN30 re-score.
- `pvt_stable_B_corrected/pvt_metrics.csv` and
  `pvt_stable_B_corrected/audit_manifest.json`: paired PVT sensitivity and
  hashes.
- Reproduction scripts are under `PAPER_NEW/scripts/` and are listed in the
  audit manifest.

Initial four-column B/C outputs made with a mis-conjugated QR row projection
are explicitly marked as superseded in
`SUPERSEDED_OUTPUTS_DO_NOT_USE.md`. No original receiver output, response
bank, observation, RINEX, solution, correction schedule, or PVT report was
overwritten.
