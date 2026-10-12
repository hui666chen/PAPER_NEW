# External-baseline status

Last updated: 2026-10-12

## Completed in this repository

- Frozen EML spacing: 0.10 chip.
- Frozen normalized Double-Delta: inner 0.10 chip, outer 0.50 chip.
- Reproducible MEDLL-like variable-projection fit.
- MEDLL ideal triangular-template versus receiver-chain-template ablation.
- Deterministic delay/amplitude/phase/C/N0 design generator.
- K-of-N gate operating-point/ROC evaluator.
- Separate algorithmic timing and end-to-end CUDA real-time-factor schema.
- One implementation smoke run on real retained PRN9 scans is summarized in
  `PILOT_SMOKE_20261012.json`; it is explicitly not a PVT result.

## Completed observation-level work

- The four external estimators execute on the retained PRN9 1110--1114 batch
  for clean, amplitude-ratio 0.05, and amplitude-ratio 0.10 scenes.
- This is 5 seeds × 3 scenes × 16 one-second windows; the compact audit has
  60 input records and 240 method schedules.
- A clean-window zero-point calibration is applied per seed and method before
  forming the correction schedules.  This removes the common prompt/tracking
  offset without using truth.
- The output is summarized in the external working directory, not committed
  as raw data:
  `gps_l1_prn09_comparable_v3/baseline_external_20261012/EXTERNAL_BASELINE_SUMMARY.json`.

A common-observable PVT replay was also completed from those schedules.  The
compact aggregate is in
`PRN09_EXTERNAL_BASELINE_PVT_20261012.json`.  Its 3-D RMSE means (metres) are:

| scene | single | EML | Double-Delta | MEDLL-ideal | MEDLL-chain |
|---|---:|---:|---:|---:|---:|
| clean | 2.047 | 2.101 | 2.128 | 3.118 | 3.850 |
| amplitude ratio 0.05 | 4.389 | 2.140 | 2.146 | 4.779 | 4.614 |
| amplitude ratio 0.10 | 7.774 | 2.347 | 2.241 | 7.135 | 9.640 |

The PVT availability means for MEDLL-chain were 0.5625, 0.6125, and 0.5875
in the three scenes; missing solutions are reported as unavailable rather than
imputed.  These are valid common-observable baseline results, but not yet the
final per-epoch estimator comparison requested for the paper: the retained
confirmation CSVs contain one 401-point window average per second.

The diagnostic supplements are summarized in
`PRN09_CRB_ROC_RUNTIME_20261012.json`:

- The 5,040-cell CRB grid gives the same local `tau`/`delta` CRB values for
  M0, M1, and M2 in the current profiled implementation; the added nuisance
  dimensions increase the condition number.  This is evidence against a
  local-Fisher-information explanation of the global jumps, not evidence that
  M1/M2 are better estimators.
- At the frozen threshold, the one-window gate has TPR 0.78125/FPR 0.0625;
  3-of-5 has TPR 0.71875/FPR 0 on one retained seed.
- The current CUDA receiver logs show 43.32--53.70 s for 145 s of signal,
  mean real-time factor 3.126.  EML and Double-Delta are sub-millisecond per
  epoch in the Python reference; MEDLL-ideal is 0.1385 s/epoch in that
  reference.  These Python values are not GPU measurements.

## Required next run

1. Retain a fresh per-epoch response batch under a persistent project path
   (not `/tmp`) for the same five seeds, or explicitly label the result as
   window-level.
2. Re-run `correlator_baselines.py` on those per-epoch files and use the clean
   zero-point calibration already frozen in the protocol.
3. Check the EML/DD sign on a clean injected control before making PVT schedules.
4. Feed every schedule through the existing common-observable PVT pipeline.
5. Complete the 560-cell grid, portability, ROC, and GPU timing studies.

The raw IQ and retained response CSVs stay outside Git.  Only compact CSV/JSON
summaries and hashes should be committed after the runs are verified.
