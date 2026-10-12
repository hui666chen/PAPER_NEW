# PRN9 paired PVT comparison: gated dual-path vs EML and Double-Delta

Date: 2026-10-12
Status: matched paired re-analysis of the retained five-seed batch; no new IQ generation or GPU run was needed.

## Design audit

The gate and classical branches use the same five noise seeds, the same source observation and navigation files per seed, the same nine-satellite geometry, and the same PVT solver/configuration. Each derived RINEX stream changes only G09 C1C. The underlying common single-path PVT output is numerically identical between the two runs after excluding input-path comments. Recomputing EML and Double-Delta from the retained 401-tap window averages reproduces every stored delay exactly (maximum absolute difference 0 chip).

Each seed contributes 16 scored seconds per scene. The statistical unit is the independent noise seed (n=5). The paired effect is `baseline 3D RMSE − gated 3D RMSE`; positive values favor the gate. Confidence intervals are percentile bootstrap intervals over the five paired seed differences (200,000 resamples, seed 2026101099).

## 3D PVT RMSE

| Scene | Gated dual-path | EML | Gate improvement vs EML (EML−gate; 95% CI) | Double-Delta | Gate improvement vs DD (DD−gate; 95% CI) |
|---|---:|---:|---:|---:|---:|
| clean | 2.047 | 2.101 | +0.053 [-0.082, +0.131] | 2.128 | +0.081 [-0.046, +0.180] |
| amp_0p05 | 3.997 | 2.140 | -1.857 [-2.145, -1.586] | 2.146 | -1.851 [-2.065, -1.678] |
| amp_0p10 | 5.133 | 2.347 | -2.786 [-2.887, -2.646] | 2.241 | -2.892 [-3.035, -2.694] |

All RMSE values are metres and are means of per-seed RMSEs. For the paired-improvement columns, a negative number means the classical baseline is better.

Calibration is not fully symmetric: the classical EML/Double-Delta branches subtract a per-seed, per-method zero point estimated as the median delay over that same seed's clean windows. The gated dual-path branch uses its independently calibrated response bank and frozen existence threshold, without subtracting a per-seed clean-window median. Thus this table compares the frozen end-to-end pipelines, but it is not a pure estimator-only comparison with identical zero-point calibration information.

## Availability and interpretation

| Scene | Gate availability | EML availability | Double-Delta availability | Gate paired record vs EML / DD (better–tie–worse) |
|---|---:|---:|---:|---:|
| clean | 1.000 | 1.000 | 1.000 | 4–0–1 / 4–0–1 |
| amp_0p05 | 1.000 | 1.000 | 1.000 | 0–0–5 / 0–0–5 |
| amp_0p10 | 1.000 | 1.000 | 0.988 | 0–0–5 / 0–0–5 |

The result is clear for this test condition: the gate is comparable to the classical correlators in clean data, but EML and Double-Delta substantially reduce 3D PVT RMSE more in both multipath scenes. At amplitude ratio 0.05, the gate is worse by about 1.86 m; at 0.10, it is worse by about 2.79–2.89 m. All five paired seeds favor each classical baseline in both multipath scenes. This does not show that the GPU software receiver itself is worse; it shows that the current PRN9 gated correction has not beaten these classical multipath controls.

The comparison is matched at the IQ-derived observation, seed, geometry, and PVT levels. It remains a window-level estimator comparison: EML/Double-Delta use a coherent one-second 401-tap average, whereas the gate uses its frozen cross-fit/persistence decision on the same underlying window data. Only five seeds and one PRN9 geometry are covered, so do not generalize beyond this setup. Because calibration information is asymmetric, treat the result as a serious operational warning and not the final estimator-only ranking; a strict method-only claim needs a follow-up with one independent clean calibration set and the same zero-point convention for every method, without retuning the frozen gate.

Machine-readable values and audit details: `PRN09_GATE_VS_CLASSICAL_PAIRED_20261012.json`.
