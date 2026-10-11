# Common-observable paired PVT results

All three methods use the same RINEX epochs, the same nine satellites,
the same navigation file, and the same RTKLIB configuration. Only PRN30
C1C differs according to the frozen correction schedule.

> Status: post-hoc infrastructure-corrected analysis of the frozen final IQ;
> this fixes the invalid separate-search PVT comparison but is not a new
> independently preregistered IQ confirmation batch.

| Scene | Method | H RMSE (m) | V RMSE (m) | 3D RMSE (m) | P95 3D (m) | Gross | Epochs | Sats |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| clean | single | 1.4591 | 1.6286 | 2.1866 | 2.9168 | 0.00% | 16 | 9 |
| clean | cross_only_refined | 1.5392 | 1.8153 | 2.3800 | 3.2685 | 0.00% | 16 | 9 |
| clean | sign_persistent_refined | 1.4591 | 1.6286 | 2.1866 | 2.9168 | 0.00% | 16 | 9 |
| amp_0p01 | single | 1.5050 | 1.8472 | 2.3827 | 3.8334 | 0.00% | 16 | 9 |
| amp_0p01 | cross_only_refined | 1.5066 | 1.8414 | 2.3792 | 3.8334 | 0.00% | 16 | 9 |
| amp_0p01 | sign_persistent_refined | 1.5050 | 1.8472 | 2.3827 | 3.8334 | 0.00% | 16 | 9 |
| amp_0p05 | single | 1.7947 | 2.1671 | 2.8137 | 4.7115 | 0.00% | 16 | 9 |
| amp_0p05 | cross_only_refined | 1.7358 | 2.0820 | 2.7106 | 4.7115 | 0.00% | 16 | 9 |
| amp_0p05 | sign_persistent_refined | 1.7947 | 2.1671 | 2.8137 | 4.7115 | 0.00% | 16 | 9 |
| amp_0p10 | single | 2.1354 | 1.8934 | 2.8540 | 3.9083 | 0.00% | 16 | 9 |
| amp_0p10 | cross_only_refined | 1.8205 | 1.5195 | 2.3713 | 3.9083 | 0.00% | 16 | 9 |
| amp_0p10 | sign_persistent_refined | 2.1217 | 1.7718 | 2.7642 | 3.9083 | 0.00% | 16 | 9 |

## Paired comparison against single-path

- clean, cross_only_refined: single minus method = -0.1934 m (-8.84%).
- clean, sign_persistent_refined: single minus method = 0.0000 m (0.00%).
- amp_0p01, cross_only_refined: single minus method = 0.0034 m (0.14%).
- amp_0p01, sign_persistent_refined: single minus method = 0.0000 m (0.00%).
- amp_0p05, cross_only_refined: single minus method = 0.1031 m (3.66%).
- amp_0p05, sign_persistent_refined: single minus method = 0.0000 m (0.00%).
- amp_0p10, cross_only_refined: single minus method = 0.4827 m (16.91%).
- amp_0p10, sign_persistent_refined: single minus method = 0.0898 m (3.15%).

## Reapplied frozen primary acceptance

- clean: PASS — non-inferiority: candidate <= single + max(0.05 m, 1%).
- amp_0p01: PASS — non-inferiority: candidate <= single + max(0.05 m, 1%).
- amp_0p05: PASS — non-inferiority: candidate <= single + max(0.05 m, 1%).
- amp_0p10: PASS — improve 3D RMSE >=0.05 m and >=1%; P95/gross not worse.
