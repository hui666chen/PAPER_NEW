# Frozen PRN9 high-elevation validation results

Status: `complete` (10/10 runs).

Target: GPS PRN9 at approximately 82.54 degrees elevation. All PVT
variants use the same nine-satellite source observations; only G09 C1C
may differ.

| Scene | Single 3D RMSE (m) | Proposed 3D RMSE (m) | Mean improvement (m) | Bootstrap 95% CI (m) | Win/tie/loss |
|---|---:|---:|---:|---:|---:|
| clean | 1.8908 | 2.0446 | -0.1538 | [-0.3574, -0.0103] | 1/3/6 |
| amp_0p05 | 2.8767 | 2.7990 | +0.0777 | [0.0095, 0.1553] | 5/4/1 |
| amp_0p10 | 4.6880 | 3.7527 | +0.9353 | [0.7256, 1.1446] | 10/0/0 |

## Frozen checks

- clean_noninferiority: FAIL
- amp_0p05_noninferiority: PASS
- amp_0p05_mean_benefit: PASS
- amp_0p10_repeatable_benefit: PASS

The ten runs vary only the preregistered SignalSim noise seed. This
experiment validates transfer from PRN30 to a new high-elevation PRN
within the same fixed date and constellation geometry; it is not a
cross-date or cross-geometry validation.
