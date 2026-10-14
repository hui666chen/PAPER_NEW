# Frozen protocol: short-delay multipath, Hatch headroom and geometry-model estimator

Date frozen: 2026-10-14 (Asia/Shanghai), before any real-data run of
`scripts/temporal_multipath_test/temporal_test.py`.

## Facts this protocol starts from (measured, results/pvt_information_test)

- Noise-free paired multipath error at 0.02-0.10 chip: EML 0.6-2.8 m
  (amplitude <= 0.25) and 2.9-5.7 m (0.50); every single-curve estimator tried
  (two-path profile, single-path fit, DD) is worse than EML.
- A single 1-s curve cannot separate a short-delay reflection from a shift of
  the direct path (bandwidth-limited resolution). The interval-PVT test
  (Stage B) confirmed that curve-only information does not beat EML.
- New information is therefore required. This protocol tests time: with a
  physical reflector, the relative phase rotates as the satellite moves.

## Physical model

Horizontal reflector at height h below the antenna:
`delay(t) = 2 h sin(el(t))`, `phase(t) = 2 pi delay(t)/lambda + phi0`,
constant amplitude ratio alpha. el(t) from the broadcast ephemeris of the
pilot scenario (BRDC 2021-06-19, receiver at 37.352721 N, 121.915773 W, 20 m),
start 10:05:30 UTC, 1800 s, the 4 highest satellites above 15 deg.

## Simulation (correlator level, no IQ)

- 1-s 401-tap curves from the calibrated PRN09 response (bank 39), direct path
  at tap 0, reflection with the model above.
- Noise: band-limited correlator noise across taps, scaled so EML 1-s noise =
  0.5 m (Stage A measured 1.54 m per 100 ms within a window).
- Response mismatch: relative residual shapes of the 32 real clean 100-ms
  windows (after a single-path fit), unscaled, switched every 60 s. This is
  conservative (it also carries 100-ms noise).
- Carrier: prompt-phase multipath error plus 2 mm noise.

## Methods

- EML: 0.1-chip power early-minus-late (identical to correlator_baselines).
- HATCH: EML + Hatch filter, window 100 s (the fair single-frequency baseline;
  longer windows are limited by ionospheric code-carrier divergence).
- MODEL: causal. Every 30 s, fit (h, alpha, phi0) on the previous 120 s of
  curves (grid h 1-30 m / alpha 0.05-0.6 / phi0 24 values, then Nelder-Mead),
  with per-epoch direct delay = EML minus the model-predicted EML bias and a
  window-constant direct amplitude. Apply the predicted bias to the next 30 s.
  A fit is used only if its relative cost gain over the no-reflection model
  exceeds the 99th percentile measured on clean (alpha = 0) runs.
- MODEL_HATCH: MODEL code and model-corrected carrier into the same Hatch filter.

Scenarios: h in {3, 5, 8, 12, 20} m, alpha in {0, 0.3, 0.5}, phi0 = 180 deg.
Primary set: h in {3, 5, 8, 12}, alpha in {0.3, 0.5}.

## Decision (RMS after first window + Hatch warm-up)

1. Headroom: primary aggregate HATCH RMS >= 0.8 m, otherwise
   **STOP_NO_HEADROOM** (Hatch is already enough; no paper here).
2. Method: primary aggregate MODEL_HATCH <= 0.6 x HATCH, at least 75% of
   primary scenarios individually <= 0.8 x HATCH, and clean (alpha = 0)
   MODEL_HATCH <= HATCH + max(5%, 0.05 m). Otherwise
   **STOP_METHOD_NOT_GOOD_ENOUGH**.
3. Both pass: **PASS_GO_TO_IQ**. The script writes a SignalSim multipath
   segment list (1-s piecewise-constant delay/phase) for the best scenario;
   next is SignalSim IQ -> GPU SDR -> the same four methods in the pseudorange
   and PVT domain, then the company simulator.

## Stated limits

The reflection is simulated with the same response the estimator uses, so
this is an upper bound on MODEL performance; real mismatch shapes reduce but
do not remove that optimism. Single horizontal reflector, constant alpha, no
ionosphere, static receiver. The decision only says whether the IQ experiment
is worth running.
