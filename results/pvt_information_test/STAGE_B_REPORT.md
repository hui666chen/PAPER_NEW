# Stage B: multi-satellite PVT Monte Carlo

Decision: **STAGE_B_STOP_INTERVAL_PVT**  
Rule: primary INTERVAL/best_comparator <= 0.85 with CI upper < 1.0, and K0 INTERVAL/DD <= 1.05  
Stage A: STAGE_A_PASS; geometry: random_pdop_le_3.5

3D RMSE (m)

| scenario | DD | EML | DD_FDE | DD_DETEX | DD_DETW | MIDPOINT | INTERVAL | INTERVAL_DET | ORACLE | best comparator | INTERVAL/best [95% CI] |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| primary_K3 | 16.18 | 9.23 | 13.23 | 18.54 | 17.56 | 10.09 | 9.93 | 15.94 | 3.66 | EML | 1.075 [1.050, 1.102] |
| clean_K0 | 8.05 | 7.37 | 10.54 | 8.37 | 8.26 | 8.01 | 8.21 | 8.02 | 3.66 | EML | 1.114 [1.091, 1.139] |
| K1 | 10.90 | 7.93 | 10.94 | 11.65 | 11.21 | 8.65 | 8.59 | 10.83 | 3.60 | EML | 1.082 [1.059, 1.109] |
| K2 | 14.43 | 8.65 | 12.33 | 16.01 | 15.36 | 9.44 | 9.36 | 14.19 | 3.60 | EML | 1.082 [1.054, 1.112] |
| K4 | 17.76 | 9.83 | 14.76 | 22.16 | 20.76 | 10.67 | 10.52 | 17.72 | 3.59 | EML | 1.070 [1.046, 1.094] |
| primary_K3_lowest_elevation | 13.21 | 8.40 | 13.02 | 17.11 | 15.07 | 9.27 | 9.69 | 13.39 | 3.57 | EML | 1.153 [1.127, 1.182] |
| K3_delay_0.02 | 9.78 | 7.68 | 10.99 | 10.33 | 10.02 | 8.49 | 8.41 | 9.64 | 3.55 | EML | 1.095 [1.074, 1.118] |
| K3_delay_0.05 | 12.75 | 8.35 | 11.61 | 14.50 | 13.96 | 9.36 | 9.40 | 12.84 | 3.63 | EML | 1.126 [1.102, 1.151] |
| K3_delay_0.10 | 14.79 | 9.00 | 12.61 | 16.59 | 15.85 | 9.71 | 9.64 | 14.55 | 3.63 | EML | 1.072 [1.046, 1.099] |
| K3_delay_0.20 | 9.77 | 8.90 | 11.48 | 14.53 | 10.20 | 8.82 | 8.62 | 9.29 | 3.58 | EML | 0.969 [0.947, 0.991] |
| primary_K3_sigma_0.5 | 16.45 | 8.69 | 14.47 | 19.05 | 17.86 | 9.45 | 9.38 | 16.31 | 1.80 | EML | 1.080 [1.053, 1.107] |
| primary_K3_sigma_2.0 | 17.09 | 11.34 | 14.10 | 20.24 | 18.77 | 11.87 | 13.22 | 17.13 | 7.15 | EML | 1.166 [1.138, 1.196] |

Semi-synthetic: per-satellite errors and intervals come from Stage A (PRN09 clean windows plus known injections); geometry and extra noise are simulated. This decides whether a hardware-in-the-loop experiment is worth running; it is not itself a PVT result.
