# Stage A: per-satellite semi-synthetic interval check

Decision: **STAGE_A_PASS**  
Rule: injected coverage >= 0.90, clean coverage >= 0.90, and every (delay>=0.05, amplitude) cell >= 0.80

- Clean coverage: 0.938; clean flag rate: 0.042
- Injected coverage: 0.906; worst cell (delay>=0.05): 0.859

| delay (chip) | amp | n | coverage | median width (m) | one-sided | sign ok (|e|>1 m) | RMS DD (m) | RMS EML (m) | RMS profile point (m) | RMS after midpoint (m) | detect | grid edge |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.02 | 0.10 | 1152 | 0.949 | 14.07 | 0.135 | 0.668 | 2.91 | 2.61 | 3.62 | 4.46 | 0.068 | 0.000 |
| 0.02 | 0.25 | 1152 | 0.911 | 14.65 | 0.140 | 0.659 | 3.19 | 2.90 | 3.67 | 4.71 | 0.071 | 0.000 |
| 0.02 | 0.50 | 1152 | 0.884 | 15.82 | 0.164 | 0.681 | 7.86 | 3.79 | 4.90 | 5.83 | 0.068 | 0.000 |
| 0.05 | 0.10 | 1152 | 0.905 | 14.65 | 0.140 | 0.646 | 3.00 | 2.72 | 3.49 | 4.51 | 0.081 | 0.000 |
| 0.05 | 0.25 | 1152 | 0.873 | 16.41 | 0.137 | 0.613 | 3.63 | 3.46 | 4.54 | 5.17 | 0.104 | 0.000 |
| 0.05 | 0.50 | 1152 | 0.859 | 19.34 | 0.214 | 0.704 | 12.14 | 5.12 | 7.02 | 7.12 | 0.202 | 0.000 |
| 0.10 | 0.10 | 1152 | 0.890 | 15.24 | 0.136 | 0.617 | 2.98 | 2.76 | 3.98 | 4.72 | 0.095 | 0.000 |
| 0.10 | 0.25 | 1152 | 0.871 | 19.63 | 0.120 | 0.592 | 4.18 | 3.70 | 5.76 | 5.73 | 0.231 | 0.000 |
| 0.10 | 0.50 | 1152 | 0.958 | 19.93 | 0.255 | 0.786 | 15.67 | 6.11 | 6.61 | 7.08 | 0.577 | 0.000 |
| 0.20 | 0.10 | 1152 | 0.873 | 17.00 | 0.135 | 0.608 | 2.91 | 2.77 | 4.44 | 5.05 | 0.178 | 0.000 |
| 0.20 | 0.25 | 1152 | 0.931 | 19.34 | 0.095 | 0.593 | 3.14 | 3.79 | 5.23 | 5.38 | 0.628 | 0.000 |
| 0.20 | 0.50 | 1152 | 0.974 | 15.24 | 0.185 | 0.687 | 7.07 | 6.27 | 4.17 | 4.42 | 0.926 | 0.000 |

## Calibration (calibration windows only)

| cell | kappa | T threshold | c_DD (chip) | c_EML (chip) |
|---|---:|---:|---:|---:|
| d005_a025 | 0.4879 | 1.851 | 0.06731 | 0.06723 |
| d005_a050 | 0.1733 | 1.267 | 0.07248 | 0.07241 |
| d010_a025 | 0.4005 | 1.784 | 0.07768 | 0.07750 |
| d010_a050 | 0.4881 | 2.23 | 0.09338 | 0.09372 |

Noise realisations within a window are epoch bootstraps of the same 100 epochs; the independent unit is the clean window (4 test windows per cell). Truth for the direct-path delay is the single-path fit on the paired clean observation.
