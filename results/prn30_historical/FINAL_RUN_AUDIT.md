# Final sign-persistence run audit

Date: 2026-09-26

## Frozen status

The protocol, estimator, receiver configs, scoring code, response banks, and
noise seed were frozen before final IQ generation.  The freeze-file SHA-256 is
`bc4f4294d5c7fdbe562ee78a4a0f5317efc8f6250f5a9a8d44685b1914771eb6`.

## Final IQ integrity

- Path: `<TMP_DIR>/prn30_persistence_iq/final_confirmation_145s_iq16.bin`
- Bytes: `14,499,900,000`
- Explicit noise seed: `2026092901`
- Clipped samples: `0`
- SHA-256: `f03095a1378e9d7b361dcc2fd25703a4e81311459bbba63c04156e9cab344531`

All four frozen response extractions contained the expected `641,601` lines.
The PRN30 state run produced 144,656 records plus its header.  Two early
loss-of-lock events were logged (each event appeared twice in the console).

## Blind decisions

| Scene | Refined single-window active | Sign-persistent active |
|---|---:|---:|
| clean | 3/16 | 0/16 |
| amplitude 0.01 | 4/16 | 0/16 |
| amplitude 0.05 | 4/16 | 0/16 |
| amplitude 0.10 | 7/16 | 2/16 |

The two emitted corrections were both in amplitude 0.10:

- window 134: `+3.1914 m`, `P_cross=0.03662`;
- window 135: `+2.8403 m`, `P_cross=0.000488`.

Thus all frozen activation-safety criteria passed.

## Frozen runtime-validity failure

The frozen protocol required every method to provide exactly 160 epochs with
exactly nine observations in every score block.  This condition failed because
independent GNSS-SDR runs decoded navigation messages at different times even
with the same nine satellites explicitly assigned.

| Method/scene | Nine-observation epochs |
|---|---:|
| single, clean | 159/160 |
| single, amplitude 0.01/0.05/0.10 | 160/160 each |
| refined single-window, clean | 0/160 (eight observations) |
| refined single-window, amplitude 0.01 | 50/160 |
| refined single-window, amplitude 0.05/0.10 | 160/160 each |
| sign-persistent, all four scenes | 160/160 each |

The frozen scorer consequently had an empty clean reference cell and stopped
instead of producing a passing result.  According to the protocol, the final
run is invalid as a complete confirmatory run.  No score window was moved and
no rerun was substituted after observing this failure.

## Diagnostic PVT result on valid nine-observation epochs

These numbers do not override the frozen runtime failure.  They show the
engineering effect where the compared native-PVT outputs met the nine-satellite
condition.

| Scene | Single 3D RMSE | Sign-persistent 3D RMSE | Change |
|---|---:|---:|---:|
| clean | 2.0762 m (159 epochs) | 2.0716 m (160 epochs) | not confirmatory because counts differ |
| amplitude 0.01 | 2.0737 m | 2.0718 m | +0.0019 m improvement; no correction emitted |
| amplitude 0.05 | 2.6162 m | 2.6158 m | +0.0004 m improvement; no correction emitted |
| amplitude 0.10 | 3.2040 m | 3.1184 m | +0.0856 m improvement (2.67%) |

For amplitude 0.10, P95 changed from `5.1012 m` to `5.0996 m`, gross-error
rate remained zero, and PVT availability was 100%.  The strong-scene
performance therefore meets the frozen improvement, P95, gross-error, and
availability criteria.  The overall frozen result remains **not passed** only
because the independent receiver runs did not share an identical observable
set.

## Required infrastructure correction

The estimator and gate should not be changed again on this batch.  The next
confirmation must compare correction schedules on one common set of tracking
observables, either by replaying one observables dump into RTKLIB PVT or by
branching one receiver run into multiple PVT paths.  This removes navigation
message acquisition time as a nuisance variable without changing the signal
model, gate, thresholds, or score windows.

