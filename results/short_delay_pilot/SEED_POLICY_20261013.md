# Short-delay pilot seed policy

Frozen before validation generation on 2026-10-13.

- Calibration seed: `2026101300`.
- Validation seeds: consecutive integers `2026101310` through `2026101314`.
- Bootstrap seed: `2026101399`.
- Validation seed count: five; it will not be increased or reduced after any
  outcome is observed.
- Mixed and counterfactual-clean IQ for one validation run use the identical
  seed. No seed is retried because its estimator result is unfavourable.
- A run may be replaced only for a recorded infrastructure failure occurring
  before estimator scoring (generation failure, clipping, failed injection
  hard gate, missing receiver output, or fewer than four common non-target
  satellites). Replacement seeds, if needed, are taken in ascending order from
  `2026101320`, and both failed and replacement run records remain visible.

