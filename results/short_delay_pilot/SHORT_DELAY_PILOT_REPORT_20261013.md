# PRN9 short-delay strong-reflection pilot: final report

Completed: 2026-10-13 (Asia/Shanghai)  
Decision: **FAIL_STOP_SUPERIORITY_CLAIM**

## Scope and audit status

This pilot tested GPS L1 C/A PRN9 at two short path separations (0.05 and
0.10 chip), two coherent reflection/direct amplitude ratios (0.25 and 0.50),
and five independent validation noise seeds. Raw tracking, coherent EML,
normalised Double-Delta (DD), and the unchanged receiver-chain-consistent
shared-ratio two-path estimator were scored on the same counterfactual
pseudorange endpoint.

Seed `2026101310` failed before pseudorange scoring because the original score
windows preceded stable G09 acquisition. The presence-only audit and the sole
timing amendment are recorded in `INFRASTRUCTURE_AMENDMENT_20261013.md`. The
five scored seeds were `2026101311`--`2026101314` and the preregistered
replacement `2026101320`. All five passed the simulator log, zero-clipping,
paired-IQ injection, tracking, and common-observable hard gates.

## Primary results

Values are equal-seed means of the seed-level pseudorange RMSE in metres.
The confidence interval is the frozen 200,000-replicate paired seed-level
bootstrap interval for `RMSE_DD - RMSE_dual`.

| Delay (chip) | Ratio | Raw RMSE (m) | EML RMSE (m) | DD RMSE (m) | Dual RMSE (m) | Dual reduction vs DD | 95% CI of DD-dual (m) | Dual wins | Pass |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---:|
| 0.05 | 0.25 | 2.499 | 3.800 | 3.909 | 4.425 | -13.18% | [-0.873, -0.221] | 0/5 | No |
| 0.05 | 0.50 | 4.004 | 3.848 | 3.968 | 3.728 | +6.05% | [-0.224, +0.685] | 3/5 | No |
| 0.10 | 0.25 | 4.999 | 2.557 | 2.611 | 4.127 | -58.03% | [-1.960, -0.749] | 0/5 | No |
| 0.10 | 0.50 | 8.579 | 7.000 | 6.273 | 8.940 | -42.51% | [-2.914, -2.501] | 0/5 | No |

The preregistered pass rule required at least one cell to show at least 30%
lower mean RMSE than DD and a paired absolute-benefit interval strictly above
zero. No cell met either combined condition. The only nominally positive cell,
0.05 chip / 0.50, improved by 6.05%, not 30%, and its interval crossed zero.

## Boundary diagnostic

Every dual estimate in all four cells was marked as a boundary solution
(160/160 scored window estimates). Both cross-fit directions selected the
minimum allowed path separation, 0.04 chip, while the recovered reflection
magnitude was at the maximum grid value, 0.50. Direct-delay estimates did not
generally hit their grid bounds. Thus the primary failure is not a missing
signal injection or a receiver acquisition failure: under the frozen model and
search grid, the two-path decomposition collapsed onto a separation/magnitude
boundary and did not deliver a stable pseudorange correction.

This diagnostic does not change the preregistered decision and must not be
used to retrofit the pilot. It explains why a larger superiority experiment
with the current estimator is not justified.

## Consequence for the paper

The result rejects the planned claim that the current two-path estimator
outperforms DD/EML for short-delay, urban-strength reflections. In accordance
with the frozen rule, no full superiority scan, gate/PVT expansion, or PRN28
extension should be launched from this estimator version.

The defensible paper direction is receiver-chain-consistent response modelling:

1. show the large response-transfer improvement already established by the
   state-consistent response experiments;
2. report this preregistered short-delay pilot as an honest operating-boundary
   experiment against DD and EML;
3. analyse the separation/magnitude boundary collapse as the current inverse
   model's limitation; and
4. treat a new short-delay parameterisation or regularised model as future
   work, requiring a newly preregistered experiment rather than tuning this
   pilot post hoc.

## Reproducibility artefacts

- Frozen protocol: `SHORT_DELAY_PILOT_PREREGISTRATION_20261013.md`
- Seed policy: `SEED_POLICY_20261013.md`
- Infrastructure amendment: `INFRASTRUCTURE_AMENDMENT_20261013.md`
- Implementation errata: `IMPLEMENTATION_ERRATA_20261013.md`
- Frozen input manifest: `FROZEN_INPUTS_20261013.json`
- Machine-readable result: `SHORT_DELAY_PILOT_RESULTS_20261013.json`
- Result SHA-256: `e849f792ce5959a0dddaeb01e5e8d19d4cb3835a2d12fa0f3f63f4d444ce5000`

