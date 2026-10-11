# Unified classical receiver-chain/PVT protocol (current writing baseline)

Status: working-paper freeze for organizing the evidence on 2026-10-11.
This document does not retroactively change historical experiments. Historical
results retain their original protocol and are labelled accordingly.

## Signal and receiver

- GPS L1 C/A;
- one selected PRN per validation batch;
- current receiver binary and a PRN/window-specific response calibration;
- 401 complex correlator taps with I/Q retained;
- response-bank GPU forward shape error recorded for every calibration bank;
- calibration data and confirmation/PVT data kept independent.

## Estimator

- receiver-chain-consistent primary and delayed responses;
- complex LS/profile estimation with shared reflection ratio;
- continuous local delay refinement only within the frozen refinement rule;
- observation-only single-window existence/risk threshold;
- no truth, injected amplitude, or true delay used in gate decisions;
- no parameter changes after confirmation data generation begins.

## Persistence gate

The current unified rule is causal **3-of-5**: at least three trusted,
same-sign correction decisions in the current five-window history are required
before the correction is applied. Otherwise the receiver falls back to the
single-path observable. The rule is not retuned separately for PRN or scene.

## PVT comparison

All three branches must be derived from one source observation stream:

1. single-path baseline;
2. always-applied dual-path correction;
3. gated dual-path correction.

Only the target PRN's C1C correction may differ. The satellite set, epoch set,
navigation file, receiver time, solver binary, and PVT configuration must be
identical. A separate receiver run for each method is not accepted as the
primary comparison because satellite acquisition timing can confound the
result.

## Scenes and metrics

Minimum scenes: clean, amplitude ratio 0.05, and amplitude ratio 0.10. A weak
0.01 scene is recommended for the safety check. Each scene should use at least
five independent noise seeds; ten or twenty are preferable for the primary
confidence interval.

Report:

- horizontal, vertical, and 3D RMSE;
- P95 3D error;
- gross-error rate, with the threshold stated;
- PVT availability;
- activation count/rate;
- paired improvement against single-path;
- seed-level confidence interval and positive-improvement fraction.

## Interpretation

- clean and weak scenes: the gated result must be non-inferior to single-path;
- strong scenes: positive improvement may be claimed only when the confidence
  interval and predeclared effect-size criterion pass;
- zero activation is a safety result, not a positive-improvement result;
- boundary solutions are reported as model/identifiability limitations, not
  silently removed.

## Version warning

Older PRN30 records use 3-of-4 in their frozen protocol files. They remain
historical evidence and cannot be reported as 3-of-5 results. A final paper
table must include a protocol-version column or exclude the old batch from the
primary pooled claim.

