# Cross-fitted basin-confidence pilot protocol

Status: frozen exploratory pilot definition. Existing endpoint outcomes have
already been inspected, so this is not a confirmatory preregistration. No rule
below may be changed after viewing the cross-fitted pilot scores.

## Information boundary

The score receives only observed per-epoch 401-tap correlations, runtime epoch
weights, the frozen response dictionary, and the candidate grid. Injection
truth is unavailable to candidate discovery, confidence calculation, and
score combination. Truth is accessed only after all scores have been produced
to calculate evaluation endpoints.

## Deterministic split

Epochs with zero-based even indices form half A and odd indices form half B.
No random split is selected. Direction A->B discovers candidate basins on A
and evaluates them on B; direction B->A reverses these roles.

## Basin discovery

For each training half, run the unchanged complete two-path profile and take
its global winner as the anchor. Traverse remaining candidates in increasing
training cost and retain at most eight competing candidates satisfying both:

1. absolute direct-delay difference from the anchor exceeds 0.05 chip;
2. Euclidean `(tau, delta_tau)` distance from every already retained
   competitor exceeds 0.05 chip.

The first condition defines a receiver-relevant harmful alternative; the
second prevents one dense grid neighbourhood from occupying the active set.

## Held-out confidence

For every fixed training-selected competitor, calculate on each held-out epoch
the weighted cost difference `competitor - anchor`. From the held-out sample
mean and covariance, use a multivariate Student-t distribution with
`df = heldout_epochs - 1` to estimate the probability that at least one
competitor has lower expected cost than the anchor. This is the directional
risk `P_A_to_B` or `P_B_to_A`.

## Combination rule

The final score is fixed as

`P_cross = max(P_A_to_B, P_B_to_A)`.

Thus both directions must support a low-risk decision. In addition,
`P_cross = 1` is forced if any of the following observable conditions holds:

- either split anchor is on the direct-delay search boundary;
- the full-data anchor is on that boundary;
- either split anchor differs from the full-data direct delay by more than
  0.05 chip;
- the two split anchors differ by more than 0.05 chip;
- a held-out covariance or active-set calculation is invalid.

No averaging, minimum-risk selection, or post-result choice among combination
rules is allowed.

## N=20 handling

The same method is tested with ten held-out epochs per direction, but uses the
small-sample Student-t calculation above. On the N=20 pilot cell, the method is
declared unsupported and fixed to single-path fallback in later experiments
if either:

- fewer than 50% of trials have split/full anchors consistent within 0.05
  chip; or
- more than 10% of all trials are both two-path-worse-than-single and assigned
  `P_cross <= 0.20`.

No alternative low-epoch score will be selected from the same pilot outcomes.

## Pilot and expansion gate

Run exactly the same first three difficult confirmatory cells used for the
naive plug-in stress test, with the same first 50 deterministic trials per
cell. Expansion to all 54 cells requires all of:

1. pooled trial-level Spearman between `P_cross` and
   `(two-path squared tau error - single-path squared tau error)` >= 0.6405;
2. Spearman is positive in every pilot cell;
3. the N=20 support rules above both pass;
4. fewer than 10% of all pilot trials are false-trust trials, defined as
   two-path worse than single-path with `P_cross <= 0.20`.

If the pilot fails, do not tune thresholds or combination rules on these data.

## Slow bootstrap role

Epoch/block bootstrap with full winner reselection is a diagnostic reference,
not a competing primary method. It is run only for trials satisfying at least
one predeclared condition: `abs(P_A_to_B-P_B_to_A)>0.25`, `0.10<=P_cross<=0.30`,
or a split/full consistency failure. It cannot replace the frozen primary
score based on which result looks better.
