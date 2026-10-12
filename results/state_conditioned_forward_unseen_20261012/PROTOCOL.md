# State-conditioned forward-response validation (unseen noise seeds)

Frozen: 2026-10-12 (Asia/Shanghai), before generating the confirmation IQ.

## Question

Does the full per-epoch state-conditioned forward model predict the measured
401-tap complex correlation response on new simulator noise seeds, beyond a
direct-path-only response?

This is a forward-model validation with the simulator's scheduled delay,
amplitude, and phase as known inputs. It is not a blind multipath-estimation
or PVT validation. No sealed truth file, estimator outcome, or holdout response
is used to choose parameters.

## Frozen batches

- PRN 9: confirmation noise seeds 2026101290--2026101294.
- PRN 28: confirmation noise seeds 2026101295--2026101299.
- Development-only seeds: PRN 9 = 2026101210; PRN 28 = 2026101230. These are
  already-seen seed values and cannot enter confirmation summaries.
- Each confirmation seed is an independent simulator noise realization. The
  calibration IQ is fixed within each PRN and generated separately from the
  validation IQ at the same absolute GNSS start time.

## Signal and response construction

- GPS L1 C/A; 25 MHz complex IQ; IF = -4.58 MHz; 401 taps from -2.00 to
  +2.00 chips in 0.01-chip steps; 100 one-millisecond epochs per 1-second
  window.
- Known reflected-path schedule for the scored interval: delay 0.8 chip,
  amplitude ratio 0.10, phase 8.833159347728436 degrees, active from relative
  second 5 through (but not including) second 21.
- The simulator is asked for a nominal 21-second file and emits 20.999 seconds
  (one millisecond less); the file sample count is checked against that
  documented generator behavior. It begins at the corresponding prior absolute scenario time
  so that scored seconds 5--20 reproduce the original PRN-specific geometry:
  PRN 9 begins 116 seconds after its original scenario epoch; PRN 28 begins
  120 seconds after its original scenario epoch. Five seconds are reserved for
  receiver acquisition/tracking before the scored interval.
- Per epoch, `P` is the response of the independent clean, isolated-PRN
  calibration IQ at that validation epoch's receiver state; `S` is generated
  by applying the simulator's linear fractional-sample delay and scheduled
  complex path coefficient to the clean calibration IQ, then correlating it
  using the same receiver state; `M` is the response from the new mixed
  validation IQ. The prediction is `P + S`, with no gain fitted to `M`.
- The PRN 28 validation JSON is reserialized as multiline JSON without changing
  its non-multipath scenario fields. This is required because the simulator
  silently ignored the prior compactly formatted `multipath` array. Every
  confirmation generation log must explicitly report one loaded segment with
  the frozen PRN, interval, delay, amplitude, and phase; otherwise that seed is
  invalid and is reported as an infrastructure failure, not silently scored.

## Primary score

For every seed, pool the 16 scored windows and 100 epochs per window, and
calculate the no-fitted-gain normalized squared errors:

`E_direct = sum |M-P|^2 / sum |M|^2`

`E_full = sum |M-(P+S)|^2 / sum |M|^2`

`benefit = 1 - E_full/E_direct`.

The inferential unit is the independent noise seed, not an epoch or tap. For
each PRN, report the mean benefit across its five confirmation seeds and the
two-sided 95% Student-t interval (`df=4`). Also report all five seed values,
per-window values, and the window-averaged complex-response NRMSE for direct
and full predictions. A fitted common complex gain is diagnostic only and
cannot replace the primary no-gain score.

## Frozen pass criteria

The forward model passes for a PRN only if both conditions hold:

1. The lower endpoint of the 95% t interval for mean seed-level benefit is
   greater than zero.
2. Mean window-averaged no-gain NRMSE of `P+S` is at most 0.10.

Report PRN 9 and PRN 28 separately. Do not pool PRNs to rescue a failed PRN.
All failures, invalid seeds, and missing windows remain in the report; no
replacement seed is allowed after any holdout response has been inspected.

## Integrity

Before confirmation generation, save hashes of this protocol, runner script,
simulator, scanner, receiver, source configs, and the corrected PRN 28 config.
Keep per-epoch observation and secondary-response CSVs, generation logs,
receiver-state CSVs, and per-seed manifests. Do not read `truth_sealed`.
