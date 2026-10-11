# Latest experiment status and valid claims

## 1. Response-model control

The controlled receiver-response comparison contains 45 cases. Replacing the
simple continuous-shift template with the receiver-chain-consistent response
reduced delay RMSE from `0.006924 chip` to `0.0001455 chip`, reduced maximum
absolute error from `0.012985 chip` to `0.000255 chip`, and increased the
within-0.01-chip hit count from `40/45` to `45/45`.

This is the strongest current method-level evidence.

## 2. PRN9 comparable validation (current 3-of-5 batch)

Five independent seeds were completed using the current receiver binary, a
PRN9-specific recalibrated response bank, unified 3-of-5 persistence, and a
common-observable PVT comparison.

| Scene | Gate activity | Mean gated 3D RMSE | Mean improvement vs single | 95% CI of improvement | Win/tie/loss |
|---|---:|---:|---:|---|---:|
| clean | 0/80 | 2.047 m | 0.000 m | [0.000, 0.000] m | 0/5/0 |
| amplitude 0.05 | 57/80 | 3.997 m | 0.392 m | [-0.045, 0.841] m | 4/0/1 |
| amplitude 0.10 | 60/80 | 5.133 m | 2.640 m | [2.405, 2.876] m | 5/0/0 |

Defensible interpretation: clean safety, non-inferiority at 0.05 under the
declared margin, and repeatable positive improvement at 0.10. The 0.05 result
is not strong evidence of statistically significant improvement.

## 3. PRN28 medium-elevation v2 validation

PRN28 was tested at `45.9325°` using current-binary calibration, a common
11-satellite observation set, and five independent seeds. The 3-of-5 gate
activated `0/80` windows in clean, 0.05, and 0.10 scenes.

The gated result exactly matched single-path in this batch. Always applying the
dual-path correction degraded 3D RMSE by `0.8738 m`, `1.9667 m`, and `0.9837 m`
in the three scenes. Therefore PRN28 proves cross-geometry safety/fallback,
not medium-elevation positive gain.

## 4. PRN4 medium-elevation historical validation

PRN4 at `49.4897°` also had zero gate activation in `240` windows. It is useful
as an independent safety case but should not be pooled with PRN28 as if it were
the same geometry or protocol version.

## 5. PRN30 historical/common-observable result

The infrastructure-corrected common-observable PRN30 result reported a gated
3D RMSE improvement of `0.0898 m` (`3.15%`) at amplitude 0.10, with no clean or
weak-scene correction. However, the associated historical records use a
3-of-4 rule and some are post-hoc infrastructure corrections. They are useful
development evidence, not interchangeable with the current 3-of-5 primary
validation.

## 6. High-elevation historical batch

An earlier PRN9 high-elevation batch had clean non-inferiority failure, while
0.05 and 0.10 passed. The later PRN9 comparable batch used a harmonized
protocol and must be kept separate. This discrepancy is a protocol/version
audit item, not evidence to be silently averaged away.

## 7. Secondary outlier-theory line

The 27-cell confirmatory outlier analysis supports retaining the Gaussian-moment
predictor as the complete semi-analytic predictor: local Wilson pass `23/27`,
global Wilson pass `26/27`, with all preregistered upgrade triggers false. The
marginal saddlepoint correction did not improve the complete joint predictor.
This is suitable for a separate theory paper, not the main receiver/PVT paper.

## Current overall conclusion

The current evidence supports:

> receiver-chain-consistent modeling plus conservative gating can improve PVT in
> identifiable strong-multipath cases and safely reject unreliable corrections.

It does not support “positive improvement at every elevation and every
reflection strength.”

