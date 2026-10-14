# Frozen amendment: profile numerical-rank audit

Date: 2026-10-13 (Asia/Shanghai)  
Scope: calibration-only implementation audit and paired coarse/dense profile comparison on the already retained 160 short-delay validation windows. No new IQ and no PVT.

## Reason for amendment

Before recalibrating κ or reading the new interval-coverage results, an implementation cross-check found that the saved coarse residual surface is exactly reproduced by the legacy profiler, but not by the mathematically equivalent affine-nuisance least-squares projection evaluated with an orthogonal QR basis. For one representative candidate, the unscaled 3-column Gram eigenvalue ratios were approximately `2.64e-13` and `5.57e-13`; the legacy `pinv(Gram, rcond=1e-12)` therefore discards both weak affine-nuisance directions. The design condition number was about `1.95e6`. This is an implementation/rank-threshold discrepancy, not evidence that the receiver data changed.

Consequently, merely replacing the delay grid while retaining the legacy pseudo-inverse would not test the intended model written in the original protocol. The amendment below isolates grid resolution under a numerically stable implementation of that stated model.

## Frozen comparison

- Primary objective for both arms: exact complex least squares after orthogonally projecting out the same constant and linear tap nuisance columns; profile the common complex gain and the frozen shared complex-ratio grid.
- Arm C (stable coarse): direct delay `[-0.40, 0.40]` chip at `0.02`-chip spacing, plus the clean-calibration zero point inserted exactly.
- Arm D (stable dense): same endpoints and zero point, at `0.002`-chip spacing.
- Both arms use the same calibrated receiver response bank, 401-tap observation, separation grid `[0.04,1.20]` at `0.02` chip, and the frozen 113-point complex-ratio grid.
- Estimate one zero point per cell as the median continuous one-path delay fit across its eight clean calibration windows, using the same stable nuisance projection.
- Calibrate each arm's κ separately, using only the 32 clean calibration windows and 100 bootstrap resamples per window (100 epochs sampled with replacement). Use NumPy `default_rng(2026101316)` and the higher empirical 90th percentile. The paired arms use the same resampled epoch indices.
- After κ is fixed, score both arms on the same existing 160 validation windows. The primary comparison is stable dense versus stable coarse. Also report the previously saved legacy coarse coverage/width as historical context, not as the affine-nuisance objective.
- Keep the previously computed existence gate unchanged. Report ungated intervals to assess profile information and gated intervals separately; gate-off remains `[0,0]`. Do not tune the gate or run PVT.

## Interpretation

This is still a post-hoc diagnostic on already inspected pilot data. It can determine whether delay-grid resolution materially changes the profile under the intended LS objective; it cannot establish independent confirmation, operational gate performance, or PVT improvement.
