# PRN28 noise-level comparability audit

Date: 2026-10-10

## Scope

This is a diagnostic comparison only. It does not modify the estimator, the
response banks, the 3-of-5 gate, or the PRN28 confirmation result. No new IQ
was generated.

## Compared quantities

The early identifiability Monte Carlo used the A-held-out PRN14 residual
covariance. It was fitted from 300 training residual epochs and checked on 200
held-out residual epochs, then used as a proper-complex Gaussian covariance
model. Its covariance trace was `8.5244236197e9` in the original correlator
scale.

For a scale check, the recorded tracking-state C/N0 and Prompt magnitude were
compared between:

- early A data used by the identifiability baseline; and
- PRN28 v2 confirmation data, windows 98, 116, and 132, pooled over five
  validation seeds.

| quantity | early A q10 | early A median | early A q90 | PRN28 q10 | PRN28 median | PRN28 q90 |
|---|---:|---:|---:|---:|---:|---:|
| C/N0 (dB-Hz) | 42.2645 | 43.7229 | 45.1981 | 41.9511 | 43.2125 | 44.5508 |
| Prompt magnitude | 2.108e6 | 2.601e6 | 3.086e6 | 2.126e6 | 2.635e6 | 3.151e6 |

The median C/N0 difference is approximately `-0.51 dB`; the Prompt-magnitude
distributions are also closely matched. There is no order-of-magnitude SNR
difference that would invalidate the early boundary as a first-order noise
reference.

## Matched-noise identifiability evidence

In the early response-consistent Monte Carlo at 100 epochs and 0.8 chip
separation:

| amplitude ratio | separation RMSE (chip) | separation outlier rate | result |
|---:|---:|---:|---|
| 0.10 | 0.11137 | 0.685 | fail |
| 0.20 | 0.04934 | 0.230 | fail |
| 0.34 | 0.01812 | 0.015 | pass under the frozen criterion |

PRN28 uses 0.05 and 0.10 amplitude ratios at 0.8 chip. Therefore the
`0.10` PRN28 condition has a directly measured, approximately noise-matched
precedent for poor separation identifiability. The `0.05` condition is weaker
than the directly tested 0.10 point, but it was not itself present in the
early grid and must be described as an extrapolation rather than a direct
measurement.

## Interpretation

The noise-level check does not support the explanation that PRN28 failed only
because it had a radically different SNR. The early Monte Carlo can therefore
be cited as supporting evidence that the 0.10/0.8-chip condition is near or
inside the known non-identifiable region.

This does not prove equality of the full PRN28 residual covariance. The raw
PRN28 IQ and per-epoch residual vectors were not retained; the remaining
open-loop CSV is a 100-epoch aggregate and includes receiver-state/model
mismatch. It must not be treated as a pure thermal-noise covariance estimate.

## Decision

Do not spend a full new 100-by-100 ideal-response Monte Carlo solely to justify
the 0.10 conclusion. The existing matched-noise baseline is sufficient for
that diagnostic claim. Keep the 0.05 statement conservative. A small
PRN28-specific ideal-response injection remains optional if a direct 0.05
claim is required by a reviewer.

The next model work, if pursued, is still justified for receiver-chain
fidelity and safe decision-making, but it should not be presented as
guaranteeing positive PRN28 PVT improvement.
