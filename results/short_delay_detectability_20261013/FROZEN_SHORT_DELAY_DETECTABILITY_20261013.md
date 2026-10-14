# Frozen protocol: short-delay detectability and profile-resolution audit

Date: 2026-10-13 (Asia/Shanghai)  
Scope: existing PRN09 calibration and short-delay pilot artifacts only; no IQ generation and no PVT.

## Question

Does the retained receiver-chain 401-tap response contain enough one-window information to distinguish a short-delay reflected component from the best single-path response, and was the earlier wide bias interval materially caused by its 0.02-chip direct-delay grid?

## Inputs and timing convention

- Use the four clean calibration response banks for windows 39, 53, 67, and 81, and their corresponding clean per-epoch response CSVs.
- Before interpreting an integration-time result, read the saved epoch timestamps and sample indices. The recorded epochs are 1-ms correlator outputs; the 100-epoch average is therefore a 100-ms average, sampled within each one-second-spaced window. It is not mislabeled as a full 1-s coherent integration.
- Estimate each bank's 1-ms residual covariance after removing a fitted one-path response, complex gain, and affine tap nuisance. Use the 8 clean calibration windows × 100 epochs available for that calibration cell, center residuals within each window, and apply a fixed shrinkage-to-diagonal estimator. Report shrinkage, eigenvalue floor, and complex pseudo-covariance ratio.
- **Calibration-only method amendment before the deterministic sweep:** the initial clean-calibration residual check showed a substantial complex pseudo-covariance ratio (about 0.4–0.5), so proper-complex-noise whitening is not assumed. The primary covariance is therefore estimated and applied in the 802-dimensional real augmented representation `[Re(taps), Im(taps)]`. This changes the covariance representation, not the test grid or thresholds; the final report will identify the amendment as calibration-diagnostic-driven and exploratory.
- Test temporal dependence using normalized residual lag covariance for lags 1–20. Compute both the IID covariance of a 100-epoch mean, C1/100, and a Bartlett-HAC sensitivity covariance. Do not claim C1/100 is valid if the lag diagnostics show material dependence.

## Deterministic response sweep

- For each calibration bank, form the noiseless response `R(0) + rho*R(delay)` from the calibrated receiver-chain response model.
- Delays: 0.02, 0.05, 0.10, 0.20 chip. Amplitude ratios: 0.10, 0.25, 0.50. Relative phases: 0°, 22.5°, …, 337.5°. Also evaluate the pilot's configured 8.83° phase as a separately labeled point.
- Fit the best single-path model by generalized least squares with continuous direct delay in [-0.40, 0.40] chip; profile its complex gain and constant/linear tap nuisances.
- In the real augmented representation, report `lambda_real = r_Rᵀ Cbar_R⁻¹ r_R` directly. Use reference dimensions 1 and 3, with false-alarm probability 1e-3 and power 0.90; report both, not an assumed single degree of freedom. These are reference thresholds, not an exact GLRT calibration for the non-regular unknown-delay multipath null. (Only under proper complex noise would `lambda_real = 2*rᴴ Cbar⁻¹ r`.)
- Estimate optimistic time-to-threshold by linear scaling from the actual 100-ms averaged observation only when the IID/HAC diagnostics support that scaling; mark it explicitly as an extrapolation to 1 s.

## Local theoretical delay range

- At every deterministic sweep point, compute the profiled Fisher-information lower bound for direct delay with separation, complex reflection ratio, complex gain, and affine tap nuisances unknown.
- Report the local 95% direct-delay interval width in metres, the effective information rank, and conditioning. Treat a singular/ill-conditioned result as unbounded, not as a narrow interval.
- The predeclared informative-range criterion is a full 95% width <= 3 m. This is a local lower bound, not a global confidence guarantee.

## Existing-data refinement gate

- Revisit the 160 short-delay validation averages only for pilot cells whose actual 8.83° phase row meets both reference detection criteria (`lambda_real` at least the df=3 threshold) and local 95% width <=3 m under the conservative covariance choice.
- For any eligible cell, refine direct delay around the existing full-grid profile candidates to 0.002 chip (or continuous local least squares), recalibrate κ only on clean calibration data, and recompute coverage/width on the existing validation data. No validation outcome may be used to select κ or the refinement rule.
- Dense-refinement rule: direct delay grid `[-0.40,0.40]` chip at 0.002-chip spacing, plus the calibration-only zero point inserted exactly; separation remains `[0.04,1.20]` at 0.02 chip; retain the already frozen ρ magnitudes/phases. Estimate one clean zero point per cell as the median continuous single-path delay over its eight clean calibration windows. Recalibrate κ as the higher empirical 90th percentile of the required κ from 100 bootstrap resamples in each of the 32 clean calibration windows (100 epochs sampled with replacement), using PRNG seed 2026101316. For each resample, required κ is the smallest tolerance that retains at least one two-path candidate at the inserted clean zero point.
- For the purpose of this diagnostic, report the dense profile interval **without the existence gate** as the primary test of whether the profile contains useful bias information; separately report the unchanged gated operational interval, assigning `[0,0]` when the existing gate is off. Do not tune or replace the gate.
- If no pilot cell is eligible, stop the refinement; report the theoretical limitation and do not retune the gate or run PVT.

## Reporting

Publish covariance diagnostics, all sweep rows, the two threshold references, local range widths, and a clear pass/stop decision. This is a post-hoc diagnostic using existing validation data, not a new confirmatory result.
