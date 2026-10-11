# PRN30 temporal-persistence confirmation audit

Date: 2026-09-26

## Status

The independent confirmation batch was generated only after the calibration
batch selected the frozen causal rule `M=4, K=3`.  The estimator, threshold
`P_cross <= 0.20`, history rule, score windows, and acceptance criteria were
not changed before the first confirmation run.

The preregistered confirmation did **not** pass all checks.

## Input integrity

- Confirmation IQ: `<TMP_DIR>/prn30_persistence_iq/confirmation_145s_iq16.bin`
- Bytes: `14,499,900,000`
- Noise seed: `2026092801`
- Clipped samples: `0`
- SHA-256: `09eb6c66151bdab7894967ed9c1663c7f557aaea4aa52fc5ec8d3ecbcf8136b7`
- PRN30 state records: `144,953` plus the CSV header
- PRN30 loss-of-lock events: `0`
- Each of the four 401-tap response files: `641,601` lines

## Frozen blind decisions

| Scene | Single-window active | Frozen 3-of-4 active |
|---|---:|---:|
| clean | 4/16 | 0/16 |
| amplitude 0.01 | 4/16 | 0/16 |
| amplitude 0.05 | 8/16 | 3/16 |
| amplitude 0.10 | 9/16 | 5/16 |

The temporal rule therefore passed its two activation-safety requirements:
zero activations in clean and no more than one activation at amplitude 0.01.

## First preregistered native-PVT run

| Scene | Single 3D RMSE | Cross-only 3D RMSE | Temporal 3D RMSE | Frozen interpretation |
|---|---:|---:|---:|---|
| clean | 5.668 m | 5.607 m | 4.242 m | PVT comparison invalidated by unequal startup satellite sets; activation count passed |
| amplitude 0.01 | 1.846 m | 2.296 m | 1.844 m | safety passed |
| amplitude 0.05 | 2.378 m | 2.484 m | 2.431 m | benefit failed; temporal was 0.053 m worse than single |
| amplitude 0.10 | 3.649 m | 3.212 m | 3.381 m | benefit passed; temporal improved 0.267 m (7.32%) |

At amplitude 0.10, cross-only improved 3D RMSE by 0.437 m while temporal
improved it by 0.267 m.  The retained benefit was therefore about 61%, below
the frozen 80% requirement.  At amplitude 0.05 neither correction method was
beneficial, so both the positive-benefit and retention requirements failed.

## Native-PVT reproducibility audit

An impossible-looking result triggered this audit: the temporal schedule had
zero correction in the clean score window but its clean PVT differed from the
single-path run.  Exact reruns showed that this was a receiver-startup nuisance,
not an effect of the zero-valued schedule.

- In the first clean-window runs, single, cross-only, and temporal used 5, 4,
  and 4 observations, respectively.
- In exact reruns they used 6, 4, and 4 observations, respectively.
- Their clean 3D RMSE consequently moved to 2.980, 13.784, and 23.229 m.
- In every later score window all methods used 9 observations.
- Between the first and second runs, later-window 3D RMSE changed by no more
  than about 0.003 m for any method.

A post-hoc fixed-satellite diagnostic assigned the same nine visible GPS
satellites (PRNs 3, 4, 7, 8, 9, 14, 16, 27, and 30) to every method.  Navigation
message decoding still completed at slightly different times in the original
42--58 s clean window.  By 60 s all methods had nine observations.  In the
still-clean 60--65 s interval, their zero-correction 3D RMSE values were:

| Method | 3D RMSE, 60--65 s |
|---|---:|
| single | 2.1083 m |
| cross-only | 2.1075 m |
| temporal | 2.1110 m |

The maximum spread is under 0.004 m.  This confirms that the apparent clean
differences came from startup satellite availability.  The preregistered clean
PVT score must not be used as evidence for or against the correction method.
The amplitude-0.01 safety result and the two strong-scene comparisons are not
affected by this startup issue because all methods used the same nine
observations in those windows and the results reproduced closely.

## Post-hoc continuous-tau diagnostic

The frozen estimator quantized direct-path delay at 0.01 chip, approximately
2.93 m.  A diagnostic retained every frozen basin, `P_cross` value, and 3-of-4
decision, but continuously refined tau within +/-0.01 chip of the frozen grid
winner.  This was performed after confirmation outcomes were known and is not
an independent confirmation result.

- All 64 optimizations converged and did not increase their profile objective.
- At amplitude 0.05, temporal 3D RMSE changed from 2.430 to 2.402 m, still
  worse than the 2.378 m single-path baseline.
- At amplitude 0.10, temporal 3D RMSE changed from 3.382 to 3.393 m, slightly
  worse than coarse temporal refinement and still better than the 3.649 m
  single-path baseline.
- Cross-only continuous refinement reached 3.107 m at amplitude 0.10, but it
  remained harmful at amplitudes 0.01 and 0.05 and is therefore not a safe
  deployable rule by itself.

The 0.05 failure is not explained solely by the 0.01-chip tau grid.  Continuous
refinement reduces overcorrection but does not create a reliable positive PVT
benefit at that strength.

## Post-hoc sign-persistence design pilot

The confirmation diagnosis suggested one additional classical consistency
check.  The current window must pass `P_cross <= 0.20`, and at least three
trusted corrections in the causal four-window history must have the same sign
as the current continuously refined correction.  There is no hysteresis,
truth input, amplitude label, or scene-dependent threshold.

This rule selected 0, 0, 0, and 3 windows in clean, amplitude-0.01,
amplitude-0.05, and amplitude-0.10 scenes, respectively.  On the post-hoc
fixed-satellite diagnostic, the results were:

| Scene | Single 3D RMSE | Sign-persistent 3D RMSE | Change |
|---|---:|---:|---:|
| amplitude 0.01 | 1.8450 m | 1.8439 m | +0.0011 m improvement |
| amplitude 0.05 | 2.3783 m | 2.3775 m | +0.0008 m improvement |
| amplitude 0.10 | 3.6485 m | 3.4297 m | +0.2189 m improvement (6.00%) |

The clean PVT numbers from this separate receiver run remain startup-confounded,
but the rule emitted no clean corrections.  The equal nine-satellite intervals
show that zero schedules agree to within a few millimetres.

This pilot is promising because it converts the amplitude-0.05 failure into a
safe fallback while retaining a useful amplitude-0.10 gain.  It was designed
after inspecting confirmation outcomes, so it is development evidence only.
It must be frozen before, and evaluated on, a new independent IQ realization.

## Defensible conclusion

The classical receiver pipeline has demonstrated three concrete results:

1. A receiver-chain-consistent two-path fit can improve native PVT in a strong
   PRN30 multipath case.  The independently confirmed gain at amplitude 0.10
   was 0.267 m (7.32%) in 3D RMSE with the frozen temporal rule; cross-only was
   stronger but unsafe in weaker cases.
2. The frozen 3-of-4 persistence rule successfully suppressed all corrections
   in clean and amplitude-0.01 windows, eliminating the weak-scene harm seen
   with single-window gating.
3. The same rule is too conservative at amplitude 0.10 and does not provide a
   positive benefit at amplitude 0.05.  It is therefore a partial result, not a
   final universally improving receiver method.

The final method candidate is now the observation-only, continuously refined,
sign-persistent rule above.  It preserves safe fallback at weak and moderate
amplitudes while admitting stable, same-sign corrections in the strong regime.
It must be frozen and tested on a new independent noise seed.  The current
confirmation batch may be used for diagnosis/design only and must not be
reused as that final confirmation.
