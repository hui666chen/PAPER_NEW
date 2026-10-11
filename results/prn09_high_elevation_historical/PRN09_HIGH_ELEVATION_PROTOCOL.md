# Frozen-design protocol for high-elevation PRN9 external validation

Design date: 2026-09-26 (Asia/Shanghai)

## Purpose and scope

This experiment evaluates the already developed classical GNSS multipath
correction on a new satellite channel without changing its estimator or
decision thresholds. The original confirmation target was PRN30 at an
elevation of approximately 27.0038 degrees. The new target is PRN9 at an
elevation of approximately 82.5378 degrees in the same fixed nine-satellite
GPS geometry.

This is an external validation across PRN and elevation, not across date,
receiver position, constellation geometry, signal type, or RF front end.

## Frozen method

The following elements are unchanged from the PRN30 V3 method:

- GPS L1 C/A, 25 MHz IQ16, IF -4.58 MHz;
- 401 complex correlator taps from -2 to +2 chips in 0.01-chip increments;
- receiver-chain-consistent response model built from an independent,
  noise-free isolated-primary calibration recording;
- complex two-component profiled least squares with constant and linear
  nuisance columns;
- the existing global dictionary and continuous local tau refinement within
  +/-0.01 chip;
- observation-only cross-fitted confidence threshold `P_cross <= 0.20`;
- causal four-window history, current-window activation, and at least three
  trusted estimates with the same correction sign (`M=4`, `K=3`);
- no hysteresis, no future-window use, and no truth-dependent decision;
- correction mapping `correction_m = -tau_hat * c / 1.023e6`;
- one common nine-satellite observation stream, with the compared RINEX files
  allowed to differ only in the target satellite's C1C field;
- the same RTKLIB PVT executable and configuration.

Only target-dependent calibration is allowed: PRN9 receives its own isolated
clean IQ, tracking-state file, and response banks. PRN30 response banks must
not be reused.

## PRN9 calibration

- Target satellite: GPS PRN9 (`G09`).
- Geometry reference epoch: 2021-06-19 10:05:40 UTC.
- Elevation: 82.537765663968 degrees.
- Azimuth: 38.291099197259 degrees.
- Calibration IQ: isolated PRN9, 145 s, noise floor -300 dB, explicit seed
  `2026100900`.
- Response-bank windows: 68, 108, and 129 s.
- The calibration IQ and response banks are generated before the validation
  freeze. They may be inspected only for processing-chain consistency and
  forward-model validation; no validation IQ or PVT outcome may exist then.
- Each response bank must pass the independent GPU zero-delay shape check with
  relative error no greater than 0.01.

## Validation scenes

Each 145 s validation realization uses the same fixed date, receiver position,
nine GPS satellites, 0.8-chip reflection delay, and 127-degree reflection
phase. SignalSim elevation power adjustment remains disabled.

| Scene | Score/extraction seconds | Multipath interval | Amplitude ratio | Response-bank window |
|---|---:|---:|---:|---:|
| clean | 61--76 | none | 0 | 68 |
| amp_0p05 | 100--115 | 98--118 | 0.05 | 108 |
| amp_0p10 | 121--136 | 119--139 | 0.10 | 129 |

Every scene contains sixteen one-second windows and 100 correlator epochs per
window. Temporal history is reset at every nonconsecutive scene boundary.

## Independent realizations and seed policy

- Planned validation runs: 10.
- SignalSim noise seeds: consecutive integers `2026100901` through
  `2026100910`, inclusive.
- Bootstrap seed: `2026100999`.
- Adaptive seed replacement, omission, or rerunning based on outcomes is
  forbidden.
- All ten valid runs are required for the confirmatory summary. Infrastructure
  failures remain associated with their original seed and are reported.

## Compared methods

1. `single`: unchanged source pseudorange.
2. `cross_only_refined`: continuously refined correction whenever the frozen
   cross-fitted confidence condition passes; retained as a secondary ablation.
3. `sign_persistent_refined`: the frozen primary PRN9 adaptation of the PRN30
   method, using the unchanged same-sign 3-of-4 rule.

The source GNSS-SDR receiver is run once per seed. All methods use the same
epochs, satellites, navigation data, and source observables. Only G09 C1C may
differ.

## Metrics and inference

For each scene and method report horizontal, vertical, and 3D RMSE; P95 3D
error; gross-error rate above 10 m; PVT availability; and activation count.
Paired improvement is `RMSE_single - RMSE_method`, so positive is better.

The independent run is the primary resampling unit. Report the mean paired
improvement over ten runs, a 100,000-resample percentile paired-bootstrap 95%
interval with seed `2026100999`, a two-sided Student-t 95% interval, standard
deviation, median, range, and win/tie/loss counts. No epoch-level bootstrap is
used for the primary claim.

Frozen interpretation criteria:

- clean is non-inferior when the lower endpoint of the primary paired-bootstrap
  interval is at least -0.05 m;
- amplitude 0.05 is non-inferior under the same rule and is called beneficial
  only if its lower endpoint is above zero;
- amplitude 0.10 supports repeatable benefit only if its lower endpoint is
  above zero, mean improvement is at least 0.05 m, and mean relative
  improvement is at least 1%;
- cross-only, P95, gross-error rate, availability, and activation counts are
  secondary and are reported regardless of direction.

## Freeze and validity

After PRN9 calibration, the protocol, exact response banks, calibration audit,
simulator configuration, estimator, orchestration code, receiver templates,
navigation assistance, binaries, seed list, and scoring code are SHA-256
sealed before any validation IQ or outcome is generated.

A run is valid only when:

- all three response CSVs contain exactly `16 * 100 * 401 + 1` lines;
- all score epochs contain exactly the frozen nine satellites including G09;
- all methods contain exactly sixteen PVT epochs per scene and identical
  timestamp sets;
- the compared RINEX streams differ only in G09 C1C;
- the native GNSS-SDR log contains a `Position at ... UTC using N
  observations` line with `N >= 4`; and
- the stricter scored PVT solutions use nine satellites at every score epoch.

No method, threshold, response-bank window, seed, scene, interval, correction
mapping, confidence-interval method, or acceptance criterion may be altered
after the validation freeze.

