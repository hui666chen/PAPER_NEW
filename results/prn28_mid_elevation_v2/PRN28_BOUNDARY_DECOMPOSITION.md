# PRN28 boundary decomposition

Purpose: post-hoc, truth-free audit of the frozen PRN28 v2 cross-fit records.
This file is diagnostic only; it does not change the gate or score PVT.

## Search bounds

- tau: [-0.40, 0.40] chip
- path separation: [0.04, 1.20] chip
- shared ratio magnitude: [0, 0.50] on the current rho grid

## Overall result

The audit contains 240 windows and 480 cross-fit H1 fits.

| boundary | count | fraction |
|---|---:|---:|
| tau lower | 0 | 0.000 |
| tau upper | 0 | 0.000 |
| separation lower (0.04 chip) | 480 | 1.000 |
| separation upper | 0 | 0.000 |
| rho magnitude upper (0.50) | 480 | 1.000 |
| rho magnitude zero | 0 | 0.000 |

Joint boundary patterns:

- `separation_lower+rho_max`: 480/480 (1.000)

## Interpretation

1. The fitted tau is not at either tau boundary.
2. Every audited fit is at the minimum path separation (0.04 chip).
3. Every audited fit is also at the current rho magnitude maximum (0.50).
4. The current implementation's `cv_no_boundary` predicate checks tau and separation, but not rho magnitude. Therefore the observed gate veto is specifically explained by the separation-lower boundary; the rho boundary is an additional unreported diagnostic warning, not the current veto condition.
5. The combination of maximum shared ratio and minimum separation is a coupled-model warning. It is not evidence that changing the gate threshold will fix PRN28. The next diagnostic is to expand the rho grid only; if the optimum remains at minimum separation, the receiver-chain response model must be audited/updated before any new gate tuning.

## Per-scene counts

| scene | direction | fits | sep lower | rho max | tau boundary |
|---|---|---:|---:|---:|---:|
| amp_0p05 | A_to_B | 80 | 80 | 80 | 0 |
| amp_0p05 | B_to_A | 80 | 80 | 80 | 0 |
| amp_0p10 | A_to_B | 80 | 80 | 80 | 0 |
| amp_0p10 | B_to_A | 80 | 80 | 80 | 0 |
| clean | A_to_B | 80 | 80 | 80 | 0 |
| clean | B_to_A | 80 | 80 | 80 | 0 |
