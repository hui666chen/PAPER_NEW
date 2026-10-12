# Recommended small-paper mainline

## Recommended title

**Receiver-Chain-Consistent Short-Delay Dual-Path Modeling and Conservative
Multipath Mitigation in a GPU GNSS Software Receiver**

An application-oriented alternative is:

**Receiver-Chain-Consistent GNSS Multipath Mitigation for Pseudorange and PVT
Using a GPU Software Receiver**

## One-sentence contribution

This work shows that a classical short-delay dual-path estimator becomes
usable in a real GNSS processing chain only after the correlator response is
calibrated to that chain and correction activation is made conservative; the
result is evaluated with identical observations at the PVT layer.

## Recommended scientific story

1. A 401-point complex correlator preserves information that ordinary E/P/L
   observations do not expose.
2. The simple “one template continuously shifted in time” assumption does not
   reproduce the actual receiver response, especially for fractional-chip
   delays.
3. A receiver-chain-consistent response bank reduces this mismatch without
   changing the estimator or global searcher.
4. A shared-ratio dual-path model estimates a direct-path delay correction.
5. A causal 3-of-5 persistence gate prevents the correction from being used in
   clean or non-identifiable windows.
6. A common-observable PVT comparison separates the correction effect from
   satellite search and navigation-message timing effects.

## Suggested paper sections

### 1. Introduction

Motivate short-delay urban-canyon multipath, response mismatch, and unsafe
over-correction. State that the paper is about making a classical model usable
in a receiver chain, not inventing a new propagation law.

### 2. GPU receiver and complex correlator observations

Describe USRP IQ input, CPU/GPU division, the 401-point complex correlator,
tracking state, and the interface to the multipath estimator. Report real-time
performance separately from the multipath contribution.

### 3. Receiver-chain-consistent dual-path model

Define the two-component complex response, nuisance complex gains, shared
reflection ratio, complex LS/profile estimation, and the calibrated response
bank. Explain why a continuous shift of one frozen template is insufficient.

### 4. Conservative correction and common-observable PVT protocol

Define the correction mapping, the observation-only threshold, causal 3-of-5
rule, and the three PVT branches: single-path, always dual-path, and gated
dual-path. All branches must use the same source observations, satellite set,
epochs, navigation data, and solver.

### 5. Results

- Response-control experiment: 45 controlled cases.
- PRN9 comparable validation: clean, 0.05, and 0.10 scenes with five seeds.
- PRN28 medium-elevation validation: safety/non-inferiority and zero activation.
- Historical PRN30 result: report only with its protocol version clearly stated.
- Optional double-delta baseline if it is available under the same observation
  protocol.

## Required external baselines before submission

The compact repository now contains a frozen protocol for EML (0.1-chip),
normalized Double-Delta, and MEDLL with ideal versus receiver-chain templates
in `protocols/BASELINE_AND_GENERALIZATION_PROTOCOL.md`.  These are required
comparators, not optional claims.  Their numerical results are pending a new
retained 401-point observation batch; no baseline result should be stated as
completed until that batch is run and the same common-observable PVT flow is
used.

### 6. Boundary analysis and limitations

Discuss zero activation, boundary solutions, identifiability, response-state
dependence, and the difference between “safe fallback” and “positive gain.”

### 7. Conclusion

Conclude that the method provides conditional PVT benefit in identifiable
strong-multipath cases and safe fallback otherwise.

## Do not combine into this paper

The semi-analytic outlier-probability theory is a possible second paper. It
should not dominate this receiver/PVT paper; otherwise the manuscript will have
two unrelated contribution stories.
