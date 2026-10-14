# Short-delay detectability diagnostic

Existing calibration and validation response data only; no IQ was generated and no PVT was run.

## Timing audit

The stored per-epoch timestamps advance by 1 ms and sample starts by 25,000 samples at 25 MHz. Thus each 100-epoch averaged 401-tap profile spans 100 ms. The nominal profile starts are one second apart; this dataset does not contain a full 1-s coherent average. Any 1-s time-to-detection below is a linear, same-phase extrapolation from 100 ms.

## Reference thresholds

- df=1: required real noncentrality 20.904.
- df=3: required real noncentrality 25.935.
- These are chi-square reference thresholds for PFA=1e-3 and power=0.90, not an exact ROC calibration for the non-regular unknown-delay multipath GLRT.

## Pilot cells at the configured 8.83° phase (HAC covariance)

| Cell | λ real | Time to df=3 threshold (s) | Local 95% delay width (m) | Eligible for dense reanalysis |
|---|---:|---:|---:|---:|
| d005_a025 | 585 | 0.004433 | 0.5074 | True |
| d005_a050 | 2565 | 0.001011 | 0.4184 | True |
| d010_a025 | 2075 | 0.00125 | 0.4492 | True |
| d010_a050 | 6486 | 0.0003998 | 0.4228 | True |

Eligible pilot rows under the frozen rule: 4.

Full phase sweep, covariance diagnostics, per-window calibration fit diagnostics, and all λ/FIM values are in `detectability_grid.csv`, `temporal_lag_diagnostics.csv`, and `SUMMARY.json`.
