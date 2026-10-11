# Instructions for another AI writing the paper

## Read in this order

1. `docs/AI_HANDOFF_EXECUTIVE_SUMMARY.md` — first read: achievements, failures, limitations, and exact writing boundaries.
2. `docs/PAPER_MAINLINE.md`
3. `docs/LATEST_EXPERIMENT_STATUS.md`
4. `docs/UNIFIED_PROTOCOL_3OF5.md`
5. `docs/CONTRIBUTION_MAP.md`
6. `docs/CLAIMS_AND_LIMITATIONS.md`
7. `results/prn09_comparable_v3/RESULTS.md`
8. `results/prn28_mid_elevation_v2/PRN28_MID_ELEVATION_V2_FINAL_REPORT.md`
9. `results/prn30_historical/COMMON_OBSERVABLE_RESULTS.md`

## Writing rules

- Treat `results/prn09_comparable_v3` and `results/prn28_mid_elevation_v2` as
  current 3-of-5 evidence.
- Treat `results/prn30_historical` and
  `results/prn09_high_elevation_historical` as versioned historical evidence.
- Never merge results with different persistence rules.
- Distinguish safety/non-inferiority from positive improvement.
- Report activation rate, seed count, confidence interval, and win/tie/loss.
- Keep the outlier-theory material as a separate optional paper line.
- Do not invent missing hardware, absolute range truth, double-delta results,
  or multi-day generalization.
- Do not turn a safe fallback result into a positive-improvement result. PRN28
  is a safety/boundary case, not a gain case.
- Always state the response-bank limitation: it is calibrated for the tested
  receiver binary and tested states, not a universal per-epoch response model.
- Always state the PRN9 0.05 confidence interval crosses zero.
- Always state that the current positive PVT evidence is conditional and is
  strongest at PRN9 amplitude 0.10.

## Desired abstract logic

Problem → receiver-chain mismatch → calibrated response model → shared-ratio
fit and 3-of-5 safety gate → common-observable PVT evaluation → conditional
benefit and safe fallback → limitations.
