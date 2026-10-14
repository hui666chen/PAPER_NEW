# Superseded audit outputs — do not use for conclusions

On 2026-10-14, a complex QR projection in the first audit implementation used
`y @ Q @ Qᴴ` instead of the correct row-vector projection `y @ Q* @ Qᵀ`.
This affected stable-QR B/C estimates and any downstream PVT generated from
those estimates.

The following initial outputs are retained for traceability only and are
superseded by files with `_corrected` in their names and the
`pvt_stable_B_corrected/` directory:

- `synthetic_noiseless.json` and `synthetic_mc_200.json` are **not affected**:
  the shared-ratio QR only projects onto real constant/linear nuisance
  columns, so the incorrect conjugation is algebraically equivalent there.
  The `_corrected` rerun is retained as a checked reproducibility copy.
- `historical_four_column_noiseless.json`: its stable B/C outputs are invalid.
- `prn30_retained_four_column_rescore.csv` and
  `prn30_retained_four_column_rescore_summary.json`: stable B/C refinements and
  dependent gate status are invalid; the historical A rows match the saved
  estimator outputs but the file as a whole should not be used.
- `pvt_stable_B/`: post-hoc PVT outputs derived from the invalid stable-B
  refinements; do not interpret as results.

The original receiver outputs, historical PVT results, and all source
observations remain unchanged. Corrected re-runs include a projection
equivalence check before scoring.
