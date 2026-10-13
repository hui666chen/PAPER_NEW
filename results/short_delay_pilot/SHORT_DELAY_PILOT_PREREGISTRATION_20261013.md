# PRN9 short-delay strong-reflection pilot preregistration

Frozen date: 2026-10-13 (Asia/Shanghai)  
Status: frozen before generating or scoring any pilot validation seed.

## Decision question

Determine whether the unchanged receiver-chain-consistent shared-ratio
two-path estimator can outperform the frozen normalized Double-Delta (DD)
correlator in a regime where classical correlators are expected to be most
challenged: short path separation and moderate-to-strong coherent reflection.
This is a small go/no-go pilot, not the final paper experiment.

No result from the historical 0.8-chip batches may be used to tune this pilot.
No PVT, elevation, persistence gate, PRN28, or old-batch zero-point analysis is
part of the primary experiment.

## Frozen signal design

- Signal: GPS L1 C/A, PRN9, static geometry beginning at 2021-06-19 10:05:30 UTC.
- Sampling/IF: 25 MHz, centre frequency 1580 MHz (IF -4.58 MHz in the receiver).
- Multipath phase: 8.833159347728436 degrees for all four cells. Phase is held
  fixed to isolate the requested delay/amplitude factors; phase generalisation
  is deferred unless this pilot passes.
- Four cells:

| Cell | Separation | Reflection/direct amplitude ratio | Score windows |
|---|---:|---:|---:|
| d005_a025 | 0.05 chip | 0.25 | 12--19 |
| d005_a050 | 0.05 chip | 0.50 | 26--33 |
| d010_a025 | 0.10 chip | 0.25 | 40--47 |
| d010_a050 | 0.10 chip | 0.50 | 54--61 |

Each multipath segment starts two seconds before and ends at least two seconds
after its eight scored windows. One 68-second mixed recording contains all four
non-overlapping segments. A second 68-second recording is generated with the
same seed and identical configuration except that the `multipath` list is
absent. The latter is the counterfactual clean reference.

## Frozen data split and seeds

- Independent isolated clean calibration seed: `2026101300`.
- Five validation seeds: `2026101310`, `2026101311`, `2026101312`,
  `2026101313`, `2026101314`.
- Bootstrap seed: `2026101399`; 200,000 paired seed-level resamples.
- The calibration seed may define response banks and zero points only. It may
  not be included in validation scoring.
- The paired clean recording of each validation seed is used only to construct
  the counterfactual scoring target and injection-quality audit. It may not be
  used to calibrate a method, choose a zero point, fit a threshold, or select a
  model.

## Injection hard gates

Before any estimator score is accepted:

1. The mixed generation log must explicitly report loading all four frozen
   time-varying multipath segments, with zero clipped samples.
2. The paired clean generation log must contain no multipath-load record and
   must report zero clipped samples.
3. For every cell, same-seed `mixed IQ - clean IQ` must have non-zero RMS in a
   frozen interior sample block. Its RMS must exceed the maximum of two frozen
   outside-segment guard blocks by at least 1000 times. Failure invalidates the
   seed; no estimator result from that seed is scored.

## Frozen observations and estimators

All methods consume the same 401-tap receiver correlator response. The tap
range is -2 to +2 chip in 0.01-chip steps. Each score is based on 100 tracking
epochs within one one-second window.

- `EML`: frozen coherent early-minus-late implementation, spacing 0.10 chip.
- `DD`: frozen normalized Double-Delta implementation, inner spacing 0.10 chip
  and outer spacing 0.50 chip.
- `dual`: the existing cross-fitted shared-complex-ratio two-path estimator,
  with its existing tau, separation, magnitude, and phase grids. In particular,
  no short-delay-specific grid refinement, new prior, gate, or optimiser change
  is allowed in the pilot.
- `raw`: no correction; reported as context, not as the primary competitor.

The two-path result for a window is the mean direct delay from the two frozen
even-to-odd and odd-to-even fits. Boundary status and recovered separation and
amplitude are reported but do not alter the correction.

## Independent zero-point definition

For every method and cell, run the unchanged method on the eight matching
windows of the independent isolated clean calibration recording. The method's
zero point is the median of its eight direct-delay estimates. The same fixed
cell-specific zero point is then subtracted from all five validation seeds.

For EML/DD, a calibration estimate is valid when its discriminator search does
not hit a boundary. For `dual`, only the direct-delay boundary is relevant to
zero-point validity because path separation is undefined in a clean signal.
At least four of eight calibration windows must be valid. Otherwise the method
is declared uncalibratable for that cell, the cell fails infrastructure QC, and
no alternative zero point is substituted.

The correction applied to the observed pseudorange error is

    correction_m = -(tau_hat - tau_zero) * c / 1.023e6.

No correction sign may be selected from validation outcomes.

## Counterfactual pseudorange-truth endpoint

Run the same receiver binary independently on each mixed/clean IQ pair. At
every scored second, form

    e_raw = (C1C_G09,mixed - C1C_G09,clean)
            - median_j(C1C_j,mixed - C1C_j,clean),  j != G09,

where the median uses all non-target satellites present in both RINEX streams
and at least four such satellites are required. The second term removes the
common receiver-time offset between the paired runs. It does not use simulator
geometric truth.

All eight frozen score seconds must have G09 and at least four common
non-target satellites in both RINEX streams. Missing any required second is an
infrastructure failure for that seed/cell; available seconds are not silently
substituted and missing values are not imputed.

For method `m`, the scored residual is `e_m = e_raw + correction_m`. For every
seed/cell/method, report RMSE, MAE, median error, 95th absolute error, gross
error rate above 10 m, valid seconds, and estimator boundary rate. The primary
unit of replication is one noise seed; the eight seconds within a seed are not
treated as independent replicates.

## Frozen go/no-go rule

For each of the four cells, calculate the five seed-level DD and dual RMSEs.
Define the paired absolute benefit as `RMSE_DD - RMSE_dual` and the paired
relative benefit as `1 - RMSE_dual / RMSE_DD`. Use the frozen bootstrap seed
and 200,000 paired seed-level resamples to obtain percentile 95% intervals.

The pilot **passes** if at least one cell satisfies both conditions:

1. the equal-seed mean dual RMSE is at least 30% below the equal-seed mean DD
   RMSE; and
2. the 95% interval for the mean paired absolute benefit lies strictly above
   zero.

All four cells are reported regardless of outcome. This pilot rule is a
go/no-go engineering criterion and does not provide multiplicity-controlled
confirmatory evidence.

- Pass: freeze a full delay/amplitude/phase/C/N0 scan, then evaluate the gate,
  hybrid architecture, PVT, and PRN28 on new data.
- Fail: stop pursuing a claim that pseudorange correction beats classical
  correlators; retain receiver-chain-consistent response modelling as the
  paper's principal technical contribution.

