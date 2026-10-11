# Audit of the frozen high-elevation PRN9 validation

Audit date: 2026-09-26 (Asia/Shanghai)

This report was written after all ten frozen validation runs completed. It
does not alter the estimator, response banks, thresholds, persistence rule,
seeds, scenes, PVT solver, confidence intervals, or acceptance criteria.

## Integrity and completion

- Target: GPS PRN9 (`G09`), elevation 82.537765663968 degrees.
- Planned/completed runs: 10/10.
- Missing or replaced seeds: none.
- Infrastructure failures: none.
- Every score epoch used the same frozen nine satellites.
- All three methods used identical solution timestamps.
- Compared RINEX streams differed only in G09 C1C.
- Native GNSS-SDR PVT success evidence was present for every run.
- Post-run verification of all frozen dependency hashes passed.
- Validation IQ and large response CSV files were deleted only after each run
  passed its hard validity checks.

Freeze SHA-256:
`7f6926c0ebe6e9baf41137ed70c40655eff7a35ed5cbd1e9fb1bbbd3cd47231c`.

## Independent PRN9 response calibration

PRN9 used dedicated isolated-primary calibration IQ and three independent
receiver-state response banks. The GPU zero-delay forward checks produced
relative shape errors of:

- window 68: `2.9819e-05`;
- window 108: `3.3041e-05`;
- window 129: `3.3876e-05`.

All are far below the frozen `0.01` limit. Consequently, the validation result
cannot reasonably be attributed to accidental reuse of a PRN30 template or a
failed PRN9 response-bank construction.

## Primary 3D-PVT endpoint

Paired improvement is `RMSE_single - RMSE_sign-persistent`; positive values
favor the proposed method.

| Scene | Single mean 3D RMSE (m) | Proposed mean 3D RMSE (m) | Mean improvement (m) | Bootstrap 95% CI (m) | Mean relative improvement | Win/tie/loss | Frozen result |
|---|---:|---:|---:|---:|---:|---:|---|
| clean | 1.8908 | 2.0446 | -0.1538 | [-0.3574, -0.0103] | -8.03% | 1/3/6 | **FAIL: not non-inferior** |
| amplitude 0.05 | 2.8767 | 2.7990 | +0.0777 | [0.0095, 0.1553] | +2.72% | 5/4/1 | **PASS: mean benefit** |
| amplitude 0.10 | 4.6880 | 3.7527 | +0.9353 | [0.7256, 1.1446] | +19.88% | 10/0/0 | **PASS: repeatable benefit** |

PVT availability remained 100% and the greater-than-10-m gross-error rate
remained zero for all three scenes and methods.

## Horizontal and vertical decomposition

The high-elevation PRN9 channel acts mainly through the vertical/clock geometry.
For the primary sign-persistent method:

| Scene | Horizontal improvement (m) | Horizontal 95% CI (m) | Vertical improvement (m) | Vertical 95% CI (m) |
|---|---:|---:|---:|---:|
| clean | -0.0174 | [-0.0445, 0.0040] | -0.1809 | [-0.4100, -0.0130] |
| amplitude 0.05 | +0.0080 | [-0.0052, 0.0207] | +0.0840 | [0.0127, 0.1657] |
| amplitude 0.10 | +0.1022 | [0.0605, 0.1406] | +0.9608 | [0.7451, 1.1740] |

Thus, both the clean-scene harm and the multipath-scene benefit are driven
primarily by the vertical component, as expected for a satellite near zenith.

## Clean-scene failure mechanism

Across 160 clean one-second windows:

- the single-window cross-fitted rule activated 64 times (40.0%);
- same-sign 3-of-4 persistence reduced this to 12 activations (7.5%);
- the 12 emitted corrections were all negative, with median approximately
  -2.50 m and range approximately -10.83 m to -0.00014 m;
- one run received a -10.83 m clean-scene correction and its 3D RMSE worsened
  by approximately 0.972 m.

The same-sign persistence rule therefore suppresses many isolated false
decisions but does not reject a systematic same-sign spurious basin. The
observation-only `P_cross` confidence can be small for these clean false
solutions, so local decision stability is not equivalent to evidence that a
reflection exists.

The unrestricted single-window ablation confirms the role of persistence: in
clean data, its mean 3D-RMSE change was -0.5660 m with a 95% interval of
[-0.9499, -0.2828] m, substantially worse than the persistent method's
-0.1538 m. Persistence helps, but it is not sufficient for safe PRN transfer.

## Scientific interpretation

The frozen estimator transfers successfully when a moderate or strong PRN9
reflection is present: amplitude 0.05 improves mean 3D RMSE, and amplitude 0.10
improves every run by an average of 0.935 m. However, the frozen PRN30-derived
confidence-and-persistence policy does not safely transfer to clean PRN9.

This result does not invalidate the completed PRN30 confirmation. It rejects a
broader claim that the current observation-only same-sign persistence rule is
automatically safe across satellite channels after response-template
recalibration alone.

The ten PRN9 runs must now be treated as a confirmation result and may not be
used both to tune a revised clean/null decision and to validate that revision.
Any new existence/null-control mechanism must be developed on separate data,
frozen, and evaluated on another independent batch.

