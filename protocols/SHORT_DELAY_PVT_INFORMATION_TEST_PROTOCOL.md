# Frozen protocol: short-delay PVT information test (semi-synthetic)

Date frozen: 2026-10-14 (Asia/Shanghai), before any real-data run of the
scripts in `scripts/pvt_information_test/`.

## Why this test

- Long delay (0.8 chip): EML and Double-Delta (DD) already bring 3D PVT RMSE
  back to the clean level (PRN09: single 7.77 m, DD 2.24 m, EML 2.35 m,
  clean 2.05 m). There is no headroom left for a new long-delay method.
- Short delay (0.05-0.10 chip): DD/EML pseudorange RMSE stays 2.6-7.0 m, so
  headroom exists, but the single-satellite two-path point estimate collapses
  to the grid boundary (160/160) and loses to DD.
- Remaining hypothesis: the 401-tap curve does not pin down one value, but it
  bounds the DD error (an interval, often one-sided), and several satellites
  together can use those bounds in the PVT.

This test checks that hypothesis cheaply, before any new simulator/IQ run.

## Stage A (per satellite, real clean data + known injected reflections)

- Data: the PRN09 short-delay pilot clean calibration windows (4 cells x 8
  windows x 100 one-millisecond epochs) and their response banks.
- Split: within each cell, even-position windows = calibration, odd-position
  windows = test. Nothing estimated on test windows is used for calibration.
- Observation: 100-epoch mean and epoch bootstraps (calibration: 24 per
  window; test: 8 per window).
- Injection: `y_mix = y_clean + g * rho * R(delta - tau1)`, where `g, tau1` are
  the single-path fit of the same clean observation. Design: delay
  {0.02, 0.05, 0.10, 0.20} chip x amplitude {0.10, 0.25, 0.50} x phase
  {0, 45, ..., 315} deg.
- Truth: direct-path tap position = single-path fit of the paired clean
  observation. DD error `e_DD = (tau1_clean - tau_DD + c_DD) * chip_m`.
- Profile: stable QR projection of constant/linear tap nuisance, direct delay
  [-0.40, 0.40] chip at 0.002 chip, separation [0.04, 1.20] at 0.02 chip, the
  existing 113-point shared-ratio grid. No raw-Gram pseudoinverse.
- Calibration (calibration windows only):
  `c_DD`, `c_EML` = median(tau_method - tau1) on clean observations;
  detection threshold = 99th percentile of `T = RSS_1path / RSS_2path - 1`
  on clean observations;
  `kappa` = 95th percentile (higher) of the smallest tolerance whose retained
  direct-delay hull contains the paired clean direct delay, over calibration
  clean observations and calibration injections (base + 5 bootstraps).
- Interval for the DD error: `[(S_min - tau_DD + c_DD), (S_max - tau_DD + c_DD)] * chip_m`.

Stage A decision: PASS only if injected coverage >= 0.90, clean coverage
>= 0.90 and every (delay >= 0.05, amplitude) cell coverage >= 0.80.
Otherwise **STOP**: the intervals are not trustworthy and Stage B is not run
(except as an explicitly labelled exploratory run).

## Stage B (multi-satellite PVT Monte Carlo)

- Geometry: random 9-satellite skies (elevation >= 10 deg, PDOP <= 3.5), or the
  real scenario sky read from an RTKLIB solution-status file.
- Error per satellite: Stage A `e_DD` (clean case for unaffected satellites,
  injected case for affected ones) plus extra noise
  `N(0, s^2 (1 + 1/sin^2 el) / 2)`, s = 1.0 m in the primary scenario.
- All satellites carry intervals; the method does not know which are affected.
- Comparators (all use the same curves): DD, EML, DD + RAIM FDE (Pfa 1e-3),
  DD + curve-detector exclusion, DD + curve-detector de-weighting (sigma x 5).
- Proposed: interval-constrained WLS (`lo_i <= b_i <= hi_i` for every
  satellite, tiny ridge toward b = 0 to make the solution unique).
- Primary scenario: N = 9, K = 3 affected, random placement, delays
  {0.05, 0.10} chip, amplitudes {0.25, 0.50}. 4000 trials, fixed seed.

Stage B decision: PASS only if
1. INTERVAL 3D RMSE <= 0.85 x the best comparator in the primary scenario,
   with the bootstrap 95% upper bound of that ratio < 1.0; and
2. in the K = 0 scenario, INTERVAL 3D RMSE <= 1.05 x DD.

Otherwise **STOP** the interval-PVT direction. Secondary scenarios (K = 1, 2,
4; lowest-elevation placement; per-delay; s = 0.5 and 2.0 m) are reported but
do not change the decision.

## What PASS / STOP mean

- PASS: worth spending simulator time on a hardware-in-the-loop experiment
  with short-delay multipath on 3-4 satellites, phases spread over 0-360 deg,
  compared against the same comparators with real RINEX/RTKLIB.
- STOP: single-frequency 1-s curves do not add PVT-useful information beyond
  DD + detection in the short-delay regime. The next candidate information
  source is multi-frequency (wideband L5/E5a delay as a constraint on L1),
  checked first by a separate theory (FIM) test.

## Limits stated in advance

Semi-synthetic: reflections are added to real clean curves, so the
injected component is an exact copy of the calibrated response (no
response mismatch for the reflection itself). Bootstrap noise within a window
is not independent; the independent unit is the clean window. Only PRN09
banks are used for every satellite in Stage B. This test can only stop or
green-light the next experiment; it is not a paper result by itself.
