# Contribution map: standard technology versus researcher-specific work

This is an academic attribution guide, not a claim that standard algorithms
were invented here.

## Excluded standard components

- GNSS signal, code, carrier, and pseudorange principles;
- GPS L1 C/A and ordinary GNSS-SDR processing stages;
- USRP X310 and CUDA as technologies;
- least squares, nonlinear least squares, profile likelihood, WLS PVT;
- the physical direct-path plus reflected-path model;
- generic temporal persistence or K-of-M decision rules;
- bootstrap, CRB, BIC, and standard statistical diagnostics.

## Researcher-specific contribution

### Core contribution — high importance

The response of the actual receiver chain is calibrated and represented as a
receiver-consistent response bank instead of assuming that a single frozen
template can simply be continuously shifted. The 45-case experiment directly
quantifies the resulting bias reduction.

### Supporting method contribution — medium importance

The calibrated 401-point complex response is connected to a shared-ratio
dual-path fit, local delay refinement, correction mapping, and a causal 3-of-5
gate in the GPU receiver workflow.

### Experimental-method contribution — medium/high importance

The common-observable PVT design isolates the effect of one PRN correction by
using the same source RINEX, epochs, satellite set, navigation data, and solver
for all branches.

### Boundary contribution — important negative evidence

PRN28 shows that zero gate activation can be the correct outcome when the
estimate is not identifiable or reaches a boundary. This defines the safe
operating region instead of claiming universal improvement.

## Conservative contribution estimate

Percentages are only a planning heuristic:

- standard theory and tools: roughly 50–60% of the technical ingredients;
- researcher-owned implementation/integration/experimental design: roughly
  25–35%;
- reviewer-recognizable method-level novelty: roughly 15–25%, concentrated in
  the response model and its safety/PVT validation.

The amount of personal engineering work is larger than the amount of formal
algorithmic novelty. That is normal for an applied receiver paper.

## What the paper should and should not claim

Should claim: a receiver-chain-consistent classical method with conditional
PVT benefit and safe fallback.

Should not claim: a new physical multipath law, a universally optimal estimator,
or a universally deployable oracle outlier probability.

