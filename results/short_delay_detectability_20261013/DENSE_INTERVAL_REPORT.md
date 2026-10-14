# Short-delay profile-grid refinement

Post-hoc diagnostic on retained response data only. No IQ was generated and no PVT was run.
Clean zero points and κ values were fixed from calibration-only data before scoring the existing validation windows.

## Numerical audit of the historical coarse profile

The saved legacy residual surface is exactly reproducible by its original profiler. However, its unscaled Gram-matrix pseudoinverse (`rcond=1e-12`) removes weak nuisance directions. The new primary comparison therefore uses stable orthogonal projection in both arms, so the coarse-vs-dense contrast isolates direct-delay grid spacing under the stated LS model.
Representative audit: `{"design_condition_number": 1948007.3102571904, "directions_below_legacy_cutoff": 2, "gram_eigenvalue_ratios_ascending": [2.6352318774922117e-13, 5.567002729525981e-13, 1.0], "legacy_pinv_rcond": 1e-12, "nuisance_column_scales": {"constant_norm": 20.024984394500787, "linear_norm": 23.180594804194893, "response_norm": 30824734.224303514}}`.
Stable affine-LS versus saved legacy surface: `{"cell": "d005_a025", "same_surface_shape": true, "saved_legacy_surface_replay_by_legacy_code": "PASS (exact equality verified independently)", "seed": 2026101311, "stable_affine_ls_vs_saved_legacy_max_relative_rss_difference": 0.8395501483179557, "stable_affine_ls_vs_saved_legacy_median_relative_rss_difference": 0.15448243791845936, "window": 36}`.

## Paired validation results

| Cell | Grid arm | κ (clean q90) | Coverage | Median width (m) | Gate active | Gated coverage | Historical legacy coverage / median width |
|---|---|---:|---:|---:|---:|---:|---:|
| d005_a025 | stable_coarse_0p02 | 0.427341 | 30/40 (75.0%) | 11.722 | 0/40 | 0/40 (0.0%) | 32/40 / 29.305 |
| d005_a025 | stable_dense_0p002 | 0.510089 | 35/40 (87.5%) | 18.169 | 0/40 | 0/40 (0.0%) | 32/40 / 29.305 |
| d005_a050 | stable_coarse_0p02 | 0.304811 | 35/40 (87.5%) | 11.722 | 0/40 | 0/40 (0.0%) | 36/40 / 29.305 |
| d005_a050 | stable_dense_0p002 | 0.369643 | 40/40 (100.0%) | 15.825 | 0/40 | 0/40 (0.0%) | 36/40 / 29.305 |
| d010_a025 | stable_coarse_0p02 | 0.248619 | 26/40 (65.0%) | 12.917 | 0/40 | 0/40 (0.0%) | 33/40 / 29.305 |
| d010_a025 | stable_dense_0p002 | 0.297101 | 38/40 (95.0%) | 17.595 | 0/40 | 0/40 (0.0%) | 33/40 / 29.305 |
| d010_a050 | stable_coarse_0p02 | 0.274929 | 34/40 (85.0%) | 12.800 | 0/40 | 0/40 (0.0%) | 37/40 / 29.305 |
| d010_a050 | stable_dense_0p002 | 0.361166 | 37/40 (92.5%) | 17.583 | 0/40 | 0/40 (0.0%) | 37/40 / 29.305 |

## Overall by grid arm

- **stable_coarse_0p02:** ungated coverage 125/160 (78.1%); median interval width 11.722 m; unchanged gate active in 0/160 windows; gated coverage 0/160 (0.0%).
- **stable_dense_0p002:** ungated coverage 150/160 (93.8%); median interval width 17.583 m; unchanged gate active in 0/160 windows; gated coverage 0/160 (0.0%).

## Interpretation

- Compare stable coarse against stable dense to assess grid spacing; do not attribute stable-vs-legacy differences to grid spacing because the legacy Gram pseudoinverse changes the effective nuisance model.
- Ungated intervals test the profile only. The unchanged gate remained off where it had been off; no PVT correction was activated by this diagnostic.
- This reuses previously inspected validation seeds/windows and is post-hoc, not independent confirmation. The 100-epoch profiles span 100 ms.
- These intervals do not establish calibrated uncertainty unless empirical coverage is adequate; neither this analysis nor the deterministic λ/Fisher calculation proves PVT improvement.
