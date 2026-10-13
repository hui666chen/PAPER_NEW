# Bias-interval PVT reanalysis rules (existing data only)

Date: 2026-10-13 (Asia/Shanghai)  
Status: rules frozen before running this reanalysis; **exploratory/post-hoc**, not a new confirmatory experiment. Earlier point-estimator results on these validation batches are already known.

## Data and scope

- Long-delay set: PRN09 comparable v3, 5 existing seeds, clean / amplitude 0.05 / amplitude 0.10, 0.8-chip reflection, one retained 401-tap average per second.
- Short-delay set: PRN09 short-delay pilot, five valid existing seeds, four cells (0.05/0.10 chip × amplitude ratio 0.25/0.50), one retained 401-tap average per scored second.
- No IQ is generated. Calibration uses only each set's independent clean calibration response data. Validation observations are not used to select κ.
- The long-delay confirmation batch has no same-seed, same-second clean RINEX counterfactual. Its residual surfaces and gate counts are analyzable, but its pseudorange-bias coverage is **not scoreable** from the retained data. No substitute truth will be fabricated.

## Grid and profile cost

- Use the receiver-chain response bank for the corresponding frozen window.
- Use the existing shared-ratio candidate construction: direct delay parameter τ in [-0.40, 0.40] chip and path separation in [0.04, 1.20] chip, both at the existing effective 0.02-chip candidate spacing; ρ=0 plus magnitudes {0.01, 0.02, 0.05, 0.10, 0.20, 0.34, 0.50} at 16 phases spaced by 22.5°.
- Save the profiled complex-LS residual for every (τ, separation, ρ) candidate for each validation window.
- A retained candidate satisfies RSS ≤ RSS_min(1+κ). The clean zero point τ₀ is estimated using the same one-path profile on the independent clean calibration response. The reported interval is the range of predicted **corrections**, b = -(τ-τ₀)×293.05 m/chip. This sign matches the existing correction interface and is compared with the required correction -e_counterfactual.

## κ calibration (fixed rule)

- For every clean calibration window, resample its existing 100 per-epoch correlator responses with replacement 100 times and average each resample to 401 taps.
- For each resample, find the smallest multiplicative tolerance κ for which the retained two-path candidates' direct-delay interval contains the calibration τ₀. If τ₀ falls between grid points, at least one retained candidate is required on each side.
- Set one κ per data family (0.8-chip family and short-delay family) to the empirical 90th percentile of these calibration-only required κ values, using the higher empirical quantile. Bootstrap PRNG seed: 2026101316.
- This is a calibration convention, not a claimed frequentist 90% coverage theorem. Validation coverage is reported as a separate test.

## Gate and scoring

- Long-delay gate: reuse the saved, truth-blind 3-of-5 gate decisions from comparable v3.
- Short-delay gate: post-hoc apply the already frozen T_null and 3-of-5 rule to the existing per-epoch responses; no validation truth enters the gate. Because the short-delay pilot had no clean validation cell, its false-positive rate is not estimable there.
- A gate-negative window is assigned [0,0] m.
- For the short-delay set, use the saved `raw_error_m` counterfactual and assess coverage of -raw_error_m. Report coverage and interval width separately by cell and delay. Use scene labels (clean vs injected) for gate false-positive/false-negative counts, not the sign or magnitude of the scored error.
- The stated ≥90% coverage criterion is applied only where an independent per-second counterfactual exists. If it fails, stop before PVT comparisons; do not retune κ on validation data.

## PVT interpretation guardrail

The requested bounded PVT comparison is not run unless step 2 passes. If reached, the bound-based solver must be validated against retained common-observation RTKLIB solutions first. For the multi-satellite semi-synthetic extension, report it as a geometry/error-propagation experiment unless every method is solved from the same underlying measurements and the absolute position truth is retained.

