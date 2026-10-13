# External-review entry point: PRN9 short-delay pilot

This directory is a self-contained audit package for a preregistered pilot
comparing an existing receiver-chain-consistent two-path estimator with
classical EML and normalised Double-Delta (DD) correlators.

## Read in this order

1. `SHORT_DELAY_PILOT_PREREGISTRATION_20261013.md` — question, estimators,
   endpoint, seeds, and pass/fail rule frozen before scoring.
2. `SEED_POLICY_20261013.md` — seed and replacement policy.
3. `INFRASTRUCTURE_AMENDMENT_20261013.md` — why the score windows were shifted
   after a presence-only acquisition audit; no pseudorange value was read.
4. `IMPLEMENTATION_ERRATA_20261013.md` — two file-size assertion corrections
   that occurred before amended validation scoring.
5. `SHORT_DELAY_PILOT_REPORT_20261013.md` — human-readable result and boundary
   diagnosis.
6. `SHORT_DELAY_PILOT_RESULTS_20261013.json` — complete machine-readable
   calibration, five seed-level results, bootstrap intervals, and decision.
7. `SHORT_DELAY_PILOT_SUMMARY_20261013.csv` — compact four-cell table.
8. `../../scripts/short_delay_pilot/run_short_delay_pilot.py` — exact runner.

## Frozen experiment

- GPS L1 C/A PRN9.
- Path separations: 0.05 and 0.10 chip.
- Reflection/direct amplitude ratios: 0.25 and 0.50.
- Fixed relative phase: 8.833159347728436 degrees.
- Five independent scored noise seeds.
- Primary endpoint: same-seed multipath-on minus multipath-off pseudorange,
  with the median common non-target-satellite difference removed.
- Primary comparison: two-path estimator versus normalised DD.
- Pass condition: at least one cell must have at least 30% lower equal-seed
  mean RMSE than DD and a paired absolute-benefit 95% bootstrap interval
  strictly above zero.

## Result that must not be reinterpreted after the fact

Decision: `FAIL_STOP_SUPERIORITY_CLAIM`.

No cell passed. The only nominally positive cell (0.05 chip, ratio 0.50)
showed 6.05% lower mean RMSE than DD, with a confidence interval crossing
zero. All 160/160 scored two-path window estimates were boundary solutions:
both cross-fit directions selected the minimum separation (0.04 chip) and the
maximum reflection magnitude (0.50).

These five validation seeds must not be reused to select a new grid,
regulariser, gate, parameterisation, or threshold. They are now test-exposed.

## Questions for an independent reviewer

Please distinguish between conclusions supported by the frozen data and ideas
that require a new experiment. In particular:

1. Does the evidence support the diagnosis of short-delay
   separation/magnitude non-identifiability or boundary collapse?
2. Could any implementation or metric choice explain the result without
   violating the frozen design?
3. Is a reparameterisation (for example, direct-delay bias plus low-order
   multipath moments), explicit model selection, or regularisation theoretically
   preferable to expanding the current `(tau, separation, rho, phase)` grid?
4. What is the smallest new, independently seeded experiment that would test
   that proposed change against DD/EML without overfitting?
5. Given the earlier strong forward-response-transfer result, is a paper centred
   on receiver-chain-consistent response modelling and empirically measured
   inverse-estimation limits technically coherent and publishable?

Please do not recommend AI/ML merely to change the paper's label; the requested
decision is about the most defensible GNSS receiver contribution supported by
the evidence.
