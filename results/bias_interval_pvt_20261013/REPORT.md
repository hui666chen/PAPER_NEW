# Bias-interval reanalysis results

Classification: exploratory/post-hoc; all inputs are retained CSV/RINEX/summary data; no IQ was generated.

- Full grid residuals retained: 400 windows × 273,347 candidates.
- κ (0.8 chip): [0.0]; κ (short delay): [0.7169136178552176].
- Comparable v3 clean gate false-positive rate: 0.000.
- Comparable v3 0.8-chip counterfactual coverage: not scoreable; no paired clean RINEX exists for the same seed and second.
- Short-delay counterfactual coverage: 0.000 (n=160 windows); median operational width 0.000 m.
- Ignoring the gate, profile-envelope coverage: 0.863; median profile width 29.305 m.
- Short-delay gate false-negative rate: 1.000.
- Stage-2 decision under the frozen ≥90% rule: **STOP_BEFORE_PVT**.

| Cell | Operational coverage | Gate active | Operational width (median m) | Ungated profile coverage | Ungated width (median m) |
|---|---:|---:|---:|---:|---:|
| d005_a025 | 0.000 | 0.000 | 0.000 | 0.800 | 29.305 |
| d005_a050 | 0.000 | 0.000 | 0.000 | 0.900 | 29.305 |
| d010_a025 | 0.000 | 0.000 | 0.000 | 0.825 | 29.305 |
| d010_a050 | 0.000 | 0.000 | 0.000 | 0.925 | 29.305 |

Since step 2 fails, the bounded-PVT and M1–M8 comparisons were not run. The full residual-surface archive (400 compressed NPZ files, approximately 774 MB) is retained in the local working directory but is not included in this GitHub upload because of its size. Per-window interval endpoints, coverage decisions, gate records, and scripts are included here. Reproducing the fits also requires the original retained correlator/RINEX input datasets, which are not part of this upload.

Files: `per_window_bias_intervals.csv`, `per_window_bias_intervals_with_profile_diagnostic.csv`, `short_delay_counterfactual_coverage.csv`, `short_delay_ungated_profile_coverage.csv`, and the five `short_gate_seed_*.json` records. The full `residual_surfaces/` archive remains local and is not uploaded.
