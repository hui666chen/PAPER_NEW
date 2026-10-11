# Frozen repeated common-observable PVT protocol, V3

Freeze date: 2026-09-26 (Asia/Shanghai)

## Reason for V3 and exclusion history

V1 seed `2026100101` was invalidated because the ninth satellite became usable
just after the relative-second-60 RINEX boundary. V2 added frozen local
broadcast-navigation assistance. This fixed ephemeris availability, and V2
seeds `2026100201` and `2026100202` completed, but seed `2026100203` showed the
remaining boundary condition: all nine ephemerides were loaded, while the G30
observable entered RINEX at relative second 61 rather than 60.

Estimator diagnostics and partial PVT outcomes from these predecessor runs
have been accessed. Consequently, V1 and V2 are infrastructure/pilot evidence
only. Seeds `2026100101`, `2026100201`, `2026100202`, and `2026100203` are
permanently excluded from the V3 confirmation claim.

V3 makes one analysis change: the clean response-extraction and score window is
moved from relative seconds 60--75 to 61--76. The frozen SignalSim multipath
schedule starts at relative second 77, so every sample in `[61, 77)` remains a
strictly clean single-path observation. The other three windows, response
banks, estimator, continuous refinement, `P_cross` threshold, causal
sign-persistence rule, geometry, multipath parameters, correction mapping, and
PVT solver are unchanged.

## Independent realizations and seed policy

- Number of planned confirmation runs: 20.
- Noise seeds: integers `2026100301` through `2026100320`, inclusive.
- Seed derivation rule: the explicit consecutive list above; adaptive seed
  selection or replacement is forbidden.
- Paired-bootstrap seed: `2026100399`.
- Reflection phase: 127 degrees for every run.
- Reflection delay: 0.8 chip.
- Amplitude ratios and score/extraction windows:
  - clean: relative seconds 61--76;
  - 0.01: relative seconds 79--94;
  - 0.05: relative seconds 100--115;
  - 0.10: relative seconds 121--136.
- Each run is a new 145 s SignalSim IQ realization using the frozen geometry.

No V3 IQ or V3 estimator/PVT outcome may be generated or inspected before this
protocol and all dependencies are hashed. No seed may be replaced, omitted, or
rerun based on its estimator or PVT outcome. Infrastructure failures remain
reported against their original seeds. Primary inference requires all 20 valid
runs; otherwise V3 is reported as incomplete.

## Frozen navigation assistance

The common-observation GNSS-SDR invocation preloads the existing files:

- `final_sign_persistence/receiver/single/gps_ephemeris.xml`;
- `final_sign_persistence/receiver/single/gps_iono.xml`;
- `final_sign_persistence/receiver/single/gps_utc_model.xml`.

Their SHA-256 hashes are included in `REPETITION_V3_FREEZE.json`. They describe
the same fixed date, geometry, and nine GPS satellites used by all runs and are
independent of V3 noise realizations.

## Frozen per-run pipeline

1. Generate the scheduled 145 s IQ realization.
2. Run the frozen PRN30 state receiver.
3. Extract sixteen 401-tap response windows for each scene, using the V3 clean
   range 61--76 and the unchanged three multipath ranges.
4. Run the unchanged observation-only estimator with the exact same window
   ranges, continuous tau refinement, `P_cross <= 0.20`, and causal same-sign
   3-of-4 rule.
5. Run GNSS-SDR once with the nine fixed satellites, frozen local navigation
   assistance, and an all-zero correction schedule to create one common RINEX
   observation stream.
6. Derive single, cross-only, and sign-persistent streams by changing only G30
   C1C according to the corresponding frozen schedule.
7. Solve all streams with the same retained RTKLIB binary, navigation file,
   and configuration.

## Hard validity checks

- Each response CSV contains exactly `16 * 100 * 401 + 1` lines.
- All score epochs contain exactly the same nine satellites.
- Every method contains exactly 16 epochs per scene.
- The three methods have identical solution timestamp sets.
- RINEX streams may differ only in G30 C1C.
- Native GNSS-SDR PVT evidence is retained; a native
  `Position at ... UTC using N observations` line with `N >= 4` establishes
  receiver PVT operation, while the stricter nine-satellite rule above governs
  V3 score validity.

The IQ and four epoch-response CSV files are hashed and deleted only after a
valid run. Compact estimator outputs, unfiltered receiver logs, RINEX, PVT
solutions, metrics, and manifests are retained.

## Metrics, confidence intervals, and interpretation

For every scene and method report horizontal, vertical, and 3D RMSE; P95 3D
error; gross-error rate above 10 m; PVT availability; and activation count.
Paired improvement is `RMSE_single - RMSE_method`, so positive is better.

The independent run is the resampling unit. Report the mean paired improvement
over 20 runs, a 100,000-resample percentile paired-bootstrap 95% interval with
seed `2026100399`, a two-sided Student-t 95% interval, standard deviation,
median, range, and positive-improvement fraction. No epoch-level bootstrap is
used for the primary claim.

- Clean, amplitude 0.01, and amplitude 0.05 are non-inferior only if the lower
  endpoint of the sign-persistent paired-bootstrap interval is at least
  `-0.05 m`.
- Amplitude 0.10 supports repeatable benefit only if that lower endpoint is
  above zero, mean improvement is at least `0.05 m`, and mean relative
  improvement is at least 1%.
- Amplitude 0.05 is called beneficial only if its lower endpoint is above zero;
  otherwise it remains a fallback/safety result.
- Cross-only, P95, gross-error, availability, and activation results are
  secondary and are always reported.

No threshold, seed, exclusion rule, confidence-interval method, window,
assistance file, or interpretation criterion may be changed after freezing.
