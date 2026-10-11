# Observation-only plug-in switch-risk protocol

Status: exploratory method-development protocol. The old/new Monte Carlo
outcomes have already been inspected in earlier analyses, so this document is
not described as a preregistration. A later independent IQ/simulator batch is
required for confirmation.

## Runtime information boundary

The score function may receive only the observed per-epoch 401-tap response,
the runtime epoch weights, the frozen response dictionary, and the frozen
noise covariance. Injected delay, amplitude, phase, true candidate index, and
true local/remote sets may be used only by the simulation wrapper and endpoint
scoring, never by the score function.

## Plug-in reference and event

The anchor is the observed global two-path profile winner. Candidate distance
is measured from this observed winner in `(tau, delta_tau)` coordinates.
Candidates within 0.10 chip are local, candidates beyond 0.20 chip are remote,
and the transition band is excluded. For computational control, the active
set consists of the eight smallest observed-cost local gaps and sixteen
smallest observed-cost remote gaps.

For each epoch, the mean response is the least-squares projection of that
epoch onto the winning model. Conditional on these fitted means, the observed
epoch weights, and the frozen covariance, exact first and second moments of
the correlated quadratic cost differences are calculated. A deterministic
Sobol Gaussian integration estimates

`P_switch = P(min remote difference < min(0, min local difference))`.

This is a plug-in stability probability, not the previously reported oracle
`P_global`.

## Mandatory sensitivity checks

Repeat selected smoke cases with active sets `(4, 8)`, `(8, 16)`, and
`(16, 32)`, and with 1,024 and 4,096 Sobol points. Do not use the score for a
deployment claim unless rankings and gate decisions are stable.

## Exploratory acceptance checks

1. Aggregate trial-level `P_switch` by its cell mean. On the new 27 cells, its
   Spearman correlation with `(two-path RMSE - single-path RMSE)` must be at
   least 0.6405, equal to 70% of the previously observed new-cell oracle
   correlation 0.915.
2. Select a gate threshold using old-27 only and apply it unchanged to new-27.
   Report the complete threshold curve and both cell-only and hierarchical
   paired bootstrap intervals.
3. Report failure cases where the selected two-path winner is wrong but
   `P_switch` is low. This directly tests the self-confirming-winner failure
   mode of plug-in inference.

Passing these checks remains exploratory because the method was designed after
viewing the existing endpoints. Confirmation requires a new frozen simulator
or IQ batch.
