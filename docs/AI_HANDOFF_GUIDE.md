# Instructions for another AI writing the paper

## Read in this order

1. `docs/PAPER_MAINLINE.md`
2. `docs/LATEST_EXPERIMENT_STATUS.md`
3. `docs/UNIFIED_PROTOCOL_3OF5.md`
4. `docs/CONTRIBUTION_MAP.md`
5. `docs/CLAIMS_AND_LIMITATIONS.md`
6. `results/prn09_comparable_v3/RESULTS.md`
7. `results/prn28_mid_elevation_v2/PRN28_MID_ELEVATION_V2_FINAL_REPORT.md`
8. `results/prn30_historical/COMMON_OBSERVABLE_RESULTS.md`

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

## Desired abstract logic

Problem → receiver-chain mismatch → calibrated response model → shared-ratio
fit and 3-of-5 safety gate → common-observable PVT evaluation → conditional
benefit and safe fallback → limitations.

