# Frozen repeated common-observable PVT protocol

Freeze date: 2026-09-26 (Asia/Shanghai)

## Purpose

Evaluate the frozen PRN30 classical multipath correction over independent
noise realizations. Each realization produces one common receiver observation
stream. The single-path, cross-only, and sign-persistent PVT solutions differ
only in the correction applied to PRN30 C1C.

This protocol contains no AI component and does not modify the estimator,
`P_cross` threshold, `M=4, K=3` sign-persistence rule, response banks, scene
windows, geometry, reflection delay, reflection phase, or PVT solver.

## Independent realizations

- Number of planned runs: 20.
- Noise seeds: integers `2026100101` through `2026100120`, inclusive.
- Reflection phase: 127 degrees for every run.
- Reflection delay: 0.8 chip.
- Amplitude ratios and score windows:
  - clean: windows 60--75;
  - 0.01: windows 79--94;
  - 0.05: windows 100--115;
  - 0.10: windows 121--136.
- Each run is a new 145 s SignalSim IQ realization using the frozen geometry.
- The previous calibration, confirmation, and final seeds are excluded.

No seed may be replaced, omitted, or rerun based on its estimator or PVT
outcome. Infrastructure failures remain reported against their original seed.
The primary analysis requires all 20 runs to pass the validity checks; an
incomplete set is reported as incomplete and is not given a primary confidence
interval.

## Frozen per-run pipeline

1. Generate one IQ file with the scheduled seed.
2. Run the frozen PRN30 state receiver and extract the four frozen 401-tap
   response blocks.
3. Run the unchanged observation-only estimator, continuous tau refinement,
   `P_cross <= 0.20`, and causal same-sign 3-of-4 rule.
4. Run GNSS-SDR once with the nine fixed satellites and an all-zero correction
   schedule to create the common RINEX observation stream.
5. Copy that stream three times and modify only PRN30 C1C according to the
   frozen zero, cross-only, or sign-persistent schedule.
6. Solve all three streams with the same retained RTKLIB `rnx2rtkp` binary,
   navigation file, and configuration.
7. Require identical solution timestamp sets and exactly 16 epochs with nine
   satellites in every scene and method.

The IQ file and the four large epoch-response CSV files are temporary. They are
hashed and deleted only after the seed has produced a valid, fully audited
result. Simulator configuration, generation log, IQ hash/size, receiver logs,
state file, compact estimator records, RINEX, PVT solutions, and manifests are
retained.

## Per-run outputs

For each scene and method report horizontal, vertical, and 3D RMSE; P95 3D
error; gross-error rate above 10 m; PVT availability; and correction activation
count. The paired improvement is

```
delta = RMSE_single - RMSE_method
```

so positive `delta` means that the correction improved PVT.

## Confidence intervals

The independent run, not the one-second window, is the resampling unit.

- Primary summary: mean paired improvement across the 20 seeds.
- Primary 95% interval: percentile paired bootstrap with 100,000 resamples and
  fixed bootstrap seed `2026100199`.
- Also report the two-sided 95% Student-t interval, sample standard deviation,
  median, minimum, maximum, and the fraction of runs with positive improvement.
- Bootstrap indices are shared across methods and metrics within a scene.
- No window-level bootstrap is used for the primary claim.

## Interpretation rules

- Clean, amplitude 0.01, and amplitude 0.05 are non-inferior only if the lower
  endpoint of the sign-persistent paired-bootstrap interval is at least
  `-0.05 m`.
- Amplitude 0.10 supports a repeatable benefit only if the lower endpoint is
  above zero, mean improvement is at least `0.05 m`, and mean relative
  improvement is at least 1%.
- Amplitude 0.05 is called beneficial only if its lower endpoint is above zero;
  otherwise it remains a fallback/safety result.
- P95, gross-error, activation, and cross-only results are secondary and are
  always reported.

No threshold, seed list, exclusion rule, confidence-interval method, or
interpretation criterion may be changed after the first new IQ is generated.
