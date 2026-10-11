# GNSS GPU Software Receiver and Receiver-Chain-Consistent Multipath Study

This is a curated, paper-oriented archive of the current GNSS software-receiver
and short-delay multipath experiments. It is intended to give another AI or
research collaborator enough controlled evidence to draft a thesis chapter,
paper outline, abstract, and figures without exposing raw IQ data or sealed
truth files.

## Current paper mainline

The recommended small-paper mainline is:

> **401-point complex correlator observations → receiver-chain-consistent
> short-delay dual-path model → observation-driven 3-of-5 conservative gate →
> common-observable PVT evaluation.**

The paper should not present GPU, NLS, or a generic K-of-M rule as new theory.
The central contribution is the receiver-specific response model and the
validated safety/benefit boundary of applying a classical dual-path correction
in a real receiver chain.

## Important status statement

The evidence supports a conditional claim:

> The calibrated dual-path correction can improve PVT in sufficiently strong
> and identifiable multipath conditions, while the temporal gate prevents
> unsafe correction when the reflection is weak, ambiguous, or the fitted
> solution reaches a boundary.

It does **not** support the claim that every satellite, elevation angle, or
reflection strength will produce a positive PVT improvement.

## Start here

1. `docs/AI_HANDOFF_EXECUTIVE_SUMMARY.md` — one-page Chinese summary of the real achievements, failures, limitations, and non-overclaim rules.
2. `docs/PAPER_MAINLINE.md` — recommended paper title, abstract logic, section outline, and figures.
3. `docs/LATEST_EXPERIMENT_STATUS.md` — latest result matrix and what each batch can prove.
4. `docs/UNIFIED_PROTOCOL_3OF5.md` — the writing/validation protocol to use going forward.
5. `docs/CONTRIBUTION_MAP.md` — standard components versus the researcher-specific contribution.
6. `docs/CLAIMS_AND_LIMITATIONS.md` — wording that is safe for a paper and wording to avoid.
7. `proposal/学硕开题报告_图表化修订版.docx` — the figure/table-enhanced opening report.

## Repository structure

- `docs/`: paper framing, contribution analysis, protocol, limitations, and AI handoff instructions.
- `results/`: compact reports and JSON summaries from the most relevant PRN batches.
- `protocols/`: reserved for final frozen protocols; historical protocols remain beside their results so that version changes are auditable.
- `scripts/`: selected analysis and experiment scripts; paths to local receiver binaries and data must be adapted before reuse.
- `figures/`: figures used in the opening report and suitable as drafts for the paper.
- `proposal/`: revised opening report.

## Data policy

Raw IQ files, large epoch-response CSV files, sealed truth, private simulator
inputs, receiver binaries, company source code, and absolute local secrets are
intentionally excluded. The compact reports retain the numerical outcomes and
protocol decisions needed for paper writing.

## Protocol warning

Some older PRN30 records use a historical 3-of-4 persistence rule, whereas the
current comparable PRN9 and PRN28 records use 3-of-5. They must not be merged
as if they were one experiment. See `docs/LATEST_EXPERIMENT_STATUS.md`.
