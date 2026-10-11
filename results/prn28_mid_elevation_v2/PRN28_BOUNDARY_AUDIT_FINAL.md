# PRN28 boundary audit — diagnostic conclusion

Date: 2026-10-10

This audit is diagnostic only. It does not modify the frozen 3-of-5 gate,
does not use sealed truth, and is not a new PVT confirmation batch.

## Frozen boundary decomposition

The existing five-seed PRN28 v2 observation contains 240 windows and 480
cross-fit fits. The fitted direct delay `tau` never touched either delay
boundary. Every fit touched both the lower separation boundary (0.04 chip) and
the old shared-ratio magnitude boundary (0.50). The current gate checks delay
and separation, not the ratio magnitude, so the formal gate veto is caused by
the separation boundary; the ratio boundary is an additional diagnostic flag.

See `PRN28_BOUNDARY_DECOMPOSITION.md` and
`PRN28_BOUNDARY_DECOMPOSITION.json`.

## Expanded-ratio diagnostic

The retained aggregate response was screened with `|rho|` extended from 0.50
to 2.0. This was not a replacement for the deleted per-epoch confirmation
files. Across 240 aggregate profiles:

| condition | count |
|---|---:|
| separation at 0.04 chip | 115 |
| separation at 1.20 chip | 28 |
| ratio at the new 2.0 limit | 13 |
| fully interior in tau, separation, and ratio | 87 |

The optimum did not move into the interior consistently. Therefore no formal
calibration or confirmation batch was started on the basis of this diagnostic.

## Receiver-chain response audit

SignalSim uses integer-sample delay followed by linear fractional-sample
interpolation. The previous `ChainResponse.exact()` implementation used a
257-tap windowed-sinc interpolation of the integer-sample response bank. On an
independent isolated-IQ forward audit, replacing the sinc operator with the
SignalSim linear operator reduced the relative response-shape error from the
percent level to approximately `4.2e-5` on average and `5.5e-5` at worst.

The corrected operator alone did not remove the PRN28 boundary. After
correcting the carrier-phase factor to use delay in samples, the retained
aggregate screening changed the fully-interior count from 87/240 to 68/240
and increased lower-separation hits from 115/240 to 157/240. This means the
interpolation bug was real, but it was not the sole cause of the PRN28
identifiability failure.

## State-conditioned response pilot

For seed `2026101130`, a clean calibration IQ was regenerated with the frozen
calibration seed. Temporary response banks were then rebuilt using the actual
PRN28 validation tracking state for representative windows 98, 116, and 132.
The calibration IQ remained independent of the mixed validation data and was
deleted after the pilot.

| scene | fixed-bank primary shape error | state-conditioned error |
|---|---:|---:|
| clean, window 98 | 5.73e-3 | 5.19e-5 |
| amplitude 0.05, window 116 | 9.29e-3 | 5.26e-5 |
| amplitude 0.10, window 132 | 3.26e-2 | 4.81e-5 |

This confirms that reusing a fixed calibration-state response bank is a real
receiver-chain mismatch. However, the corresponding expanded-ratio profiles
still reached a boundary:

| scene | state-conditioned diagnostic optimum |
|---|---|
| clean | separation 0.04 chip, `|rho|=0.50` |
| amplitude 0.05 | separation 0.04 chip, `|rho|=0.50` |
| amplitude 0.10 | separation 1.20 chip, `|rho|=0.02` |

The first state-conditioned profile printout used an incorrect scalar carrier
phase in the temporary linear operator and was discarded. The values above
come from the corrected re-score, with
`exp(j*omega_per_sample*delay_chip/chip_per_sample)`.

This pilot therefore identifies the next model requirement: the formal
estimator must use per-epoch state-conditioned response columns, rather than a
single averaged state-conditioned bank. It is not evidence that a new gate
threshold will solve the problem.

## Decision

Do not start the requested independent mid-elevation confirmation batch yet.
The required acceptance conditions — clean non-degradation and stable gate
activation with PVT improvement at both 0.05 and 0.10 — cannot be tested fairly
until the per-epoch state-conditioned model has passed a forward validation.

The next implementation is:

1. generate per-epoch clean response columns from the independent calibration
   IQ and the observed receiver state;
2. replace the common dictionary in the profile by the corresponding
   per-epoch dictionary while keeping the same delay grid, shared-ratio grid,
   cross-fit split, 3-of-5 rule, and PVT code;
3. validate the forward response on held-out clean calibration intervals;
4. only then create a fresh independent calibration batch and a fresh
   confirmation batch.

No gate, threshold, or PVT scoring rule was changed in this audit.

Corrected diagnostic outputs:

- `linear_delay_response_probe/summary.json`
- `linear_model_profile_2p0_correct_phase/summary.json`
- `state_conditioned_seed_2026101130/reprofile_correct_phase/summary.json`
