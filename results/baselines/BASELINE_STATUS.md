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

## Not yet completed

No numerical EML, Double-Delta, MEDLL, CRB-extended-grid, portability, ROC, or
GPU timing result is claimed here.  The old PRN9 confirmation runner deleted
the mixed 401-point response CSVs after use.  A fresh retained observation
batch is therefore required before the external baselines can be scored.

## Required next run

1. Run `scripts/baselines/run_prn09_baseline_observation.py` for one pilot seed.
2. Run `scripts/baselines/correlator_baselines.py` for clean, 0.05, and 0.10.
   Use `--max-windows 1` first to verify the live chain and the schedule
   adapter; remove it only for the full 16-window baseline.
3. Check the EML/DD sign on a clean injected control before making PVT schedules.
4. Process all five frozen seeds only after the pilot passes.
5. Feed every schedule through the existing common-observable PVT pipeline.

The raw IQ and retained response CSVs stay outside Git.  Only compact CSV/JSON
summaries and hashes should be committed after the runs are verified.
