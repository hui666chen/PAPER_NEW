# Implementation erratum before validation generation

Date: 2026-10-13 (Asia/Shanghai)

The first invocation stopped immediately after generating the independent
calibration IQ and before generating any validation mixed/clean IQ, running an
estimator, or accessing an outcome. The runner expected `6,799,900,000` bytes
by extrapolating an older 145-second file-size convention. SignalSim produced
exactly `6,800,000,000` bytes for the frozen 68-second configuration, equal to
`68 s × 25,000,000 complex samples/s × 4 bytes/complex sample`.

The implementation was amended only to:

1. set `EXPECTED_IQ_BYTES` to `6,800,000,000`; and
2. resume from that already generated calibration IQ when its generation log
   records zero clipping, instead of regenerating the identical file.

No signal parameter, seed, estimator, zero-point rule, endpoint, confidence
interval, or pass/fail threshold changed. The failed invocation and original
generation log remain in the raw experiment directory. The amended script is
hashed in the frozen-input manifest before validation generation begins.

## Second file-length assertion correction after the acquisition amendment

After the presence-only acquisition audit moved the scenes into a 96-second
recording, the first amended invocation again stopped at the calibration file
length assertion, before running the receiver, any estimator, or any amended
validation seed. SignalSim reported completing all 96,000 ms with zero clipped
samples and wrote `9,599,900,000` bytes (95.999 s of complex IQ16 at 25 MHz),
one millisecond below the nominal configuration duration.

The runner's exact expected size was changed only from `9,600,000,000` to the
deterministic size actually recorded in the generation log,
`9,599,900,000`. All scheduled score and guard blocks end before 94 seconds,
so the missing final millisecond cannot overlap any analysed sample. The
already generated calibration IQ is reused; no seed or signal is regenerated.
No experimental or statistical rule changed.
