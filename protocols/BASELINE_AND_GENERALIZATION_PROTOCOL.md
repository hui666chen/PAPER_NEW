# Frozen baseline and generalization protocol

Status: design frozen; numerical baseline results are pending a new retained
401-point observation batch.  The old PRN9 confirmation runner deleted the
mixed response CSVs after PVT, so no baseline number is back-filled from a
missing file.

## 1. Purpose

The current three-way comparison (single-path, always dual-path, gated
dual-path) contains only variants of the proposed processing chain.  This
protocol adds established correlator and MEDLL comparators while keeping the
same 401-point complex observations, receiver states, satellite set, epochs,
navigation data, and common-observable PVT solver.

The receiver-chain contribution is isolated by one controlled ablation:

```text
MEDLL + ideal triangular template
            versus
MEDLL + calibrated receiver-chain response bank
```

Only the template changes.  The fit range, nuisance terms, coarse grid,
refinement, stopping rules, and PVT conversion remain identical.

## 2. Baseline methods

### 2.1 EML / narrow early-minus-late

Use coherent prompt power from the 401-point scan and a fixed 0.10-chip
early-minus-late spacing.  On a delay grid `u`,

```text
D_EML(u) = |R(u-0.05)|^2 - |R(u+0.05)|^2 .
```

The reported delay is the zero crossing nearest zero; if no crossing exists,
the minimum absolute discriminator value is reported and `boundary=true` is
recorded.  The grid is `[-0.50, 0.50]` chip with 0.002-chip spacing.

### 2.2 Normalized Double-Delta

The first comparator is the BPSK normalized double-delta form, not HRC/Strobe;
HRC/Strobe is not silently mixed into the same result table.  Inner spacing is
0.10 chip and outer spacing is 0.50 chip:

```text
D_DD(u) = [P(u-0.05)-P(u+0.05)]/0.10
          - [P(u-0.25)-P(u+0.25)]/0.50 .
```

The nearest zero-crossing rule and search grid are the same as EML.  These
spacings are frozen before the baseline run; they are not tuned to a result.

### 2.3 MEDLL-like two-component fit

For each epoch, fit

```text
y(x) = a0 h(x; tau) + a1 h(x; tau-delta) + b0 + b1 x + n(x),
```

where `a0`, `a1`, `b0`, and `b1` are complex linear nuisance coefficients and
`tau`, `delta` are searched nonlinearly.  The direct-path correction uses the
same sign convention as the existing receiver estimator:

```text
correction_m = -tau_chip * 299792458 / 1.023e6 .
```

Frozen fit settings:

| setting | value |
|---|---:|
| fit interval | `[-1.50, 1.50]` chip |
| direct delay grid | `[-0.50, 0.50]`, step `0.02` chip |
| excess-delay grid | `[0.02, 1.00]`, step `0.04` chip |
| refinement starts | best 8 separated coarse candidates |
| linear terms | complex `a0`, `a1`, constant `b0`, slope `b1` |
| optimizer | bounded L-BFGS-B, maximum 400 iterations |

The ideal template is a unit triangular C/A lobe, aligned to the peak of the
zero-delay chain response for that bank.  The chain template is the calibrated
bank loaded from `bank_<window>.npz` using the bank's recorded processing
specification.  No test-scene truth is used to choose either template.

MEDLL is a classical comparator, not a claimed new algorithm.  The paper may
claim only the controlled effect of the receiver-chain template and the
additional observation-driven gate.

## 3. First retained-observation pilot

Before the large parameter sweep, regenerate a new PRN9 observation batch with
the current receiver binary and retain the mixed response CSVs:

```bash
python3 scripts/baselines/run_prn09_baseline_observation.py \
  --output-root /tmp/prn09_baseline_observations \
  --temp-root /tmp/prn09_baseline_tmp \
  --seeds 2026101110
```

The output root is intentionally outside Git.  It contains large response
CSV files and is not uploaded.

For every scene and seed, run:

```bash
python3 scripts/baselines/correlator_baselines.py \
  --scan-csv /tmp/prn09_baseline_observations/seed_2026101110/open_loop/clean/prn09_epoch_responses.csv \
  --chain-bank-dir /home/zoxoai/Documents/ChatGPT/paper/gps_l1_prn09_comparable_v3/response_banks \
  --bank-window 68 \
  --config scripts/baselines/baseline_config.json \
  --output-csv /tmp/prn09_baseline_observations/seed_2026101110/baselines_clean.csv
```

Repeat the command for `amp_0p05` and `amp_0p10`, then for the remaining four
frozen seeds.  The output is observation-only; PVT scoring is a separate stage.

## 4. PVT comparison and metrics

Each baseline must produce a correction schedule and enter the existing
common-observable PVT pipeline.  The three proposed branches remain:

1. single-path;
2. always shared-ratio dual-path;
3. current observation-driven 3-of-5 gated dual-path.

The added branches are EML, normalized Double-Delta, MEDLL-ideal, and
MEDLL-chain.  Every branch uses the same source observation stream and must
have an explicit schedule hash.  Report horizontal, vertical, 3D RMSE, P95
position error, gross-error rate, PVT availability, activation rate, and
paired win/tie/loss.  Do not pool historical 3-of-4 PRN30 results with the
current 3-of-5 batch.

The first required table is:

```text
single | always-dual | gated-dual | EML-0.10 | Double-Delta
       | MEDLL-ideal | MEDLL-chain
```

The second required table is the direct template ablation:

```text
MEDLL-ideal - MEDLL-chain
```

with identical observations and paired per-seed confidence intervals.

For the classical baselines, `scripts/baselines/baseline_schedule.py` uses a
frozen deterministic adapter: take the median valid epoch delay in each
one-second window; if fewer than 50% of that window's epoch fits are finite and
interior, write zero correction and count the window as rejected.  This
validity rule is reported separately and is not allowed to be tuned from PVT
outcomes.

## 5. Parameter coverage

The confirmation design is a Cartesian grid, generated by
`scripts/baselines/design_parameter_sweep.py`:

```text
delay: 0.05, 0.10, 0.20, 0.40, 0.60, 0.80, 1.00 chip
amplitude ratio: 0.10, 0.20, 0.34, 0.50, 0.70
relative phase: 0, 90, 180, 270 degree
C/N0: 35, 40, 45, 50 dB-Hz
```

This is 560 cells.  The pilot must pass implementation checks before the full
grid is run.  The sweep reports delay-error RMSE, 95-percentile error,
outlier rate, boundary rate, and (when a PVT geometry is supplied) horizontal,
vertical, and 3D PVT error.  The 0.80-chip case is retained as a stress case,
not labelled a typical urban short-delay case.

## 6. CRB and identifiability boundary

Reuse the verified profiled Fisher implementation in
`scripts/verify_crb_implementation.py` and the existing M0/M1/M2 evaluator.
Extend the grid with C/N0 by scaling the independently estimated 401-tap
noise covariance from its recorded reference C/N0.  For every cell report:

```text
tau CRB SD, delta-tau CRB SD, Fisher rank, condition number,
phase sensitivity, and whether the fitted search reaches a boundary.
```

Do not convert a CRB threshold into a universal physical limit.  If a binary
plot is needed for this experiment, use the predeclared diagnostic line
`delta-tau CRB SD <= 0.02 chip` and show the raw continuous surface beside it.
PRN28's zero activation is a correct rejection only if its measured state and
noise level are placed on this map and the response-control audit is passed.

## 7. Response-bank portability

The portability test is a controlled matrix, not a claim that one bank works
everywhere:

```text
PRN: 9, 28, one additional stable PRN if available
Doppler: low, nominal, high state from tracking CSV
bandwidth/sample-rate: the installed 25-MHz chain and each available
                         independently configured front-end bandwidth
```

For each cell, build the bank from independent clean calibration data, verify
the forward response, then evaluate a held-out clean/mixed response.  Report
the maximum normalized shape error and the delay bias.  A bank may not be
reused across PRNs or receiver binaries without a metadata/hash match.

For PRN28, first extend the diagnostic reflection-ratio grid beyond the old
upper bound.  If clean, 0.05, and 0.10 all land on the same ratio/separation
boundaries, label the outcome as an identifiability/model-boundary case; do
not change the gate merely to force activation.

## 8. Gate detection performance

Run `scripts/baselines/gate_roc.py` on observation-only records.  The primary
score is the existing `cv_gain_min` with the existing no-boundary condition.
Evaluate every `K` for `N in {1, 3, 5, 7}`; `K=3,N=5` remains the frozen
primary rule.  Use clean/weak scenes as negative or positive labels only in
the final scoring table, never in the gate decision.

Report single-window TPR/FPR, K-of-N TPR/FPR, delayed activation, and the
per-seed bootstrap interval.  An SQM comparison is marked **not available**
unless an independently computed SQM column is supplied from the same
observations; it must not be fabricated from `cv_gain_min`.

## 9. Computation cost

Record separately:

1. CUDA 401-tap correlation time per channel per epoch;
2. NLS/shared-ratio time per epoch;
3. gate time per epoch;
4. end-to-end wall time and real-time factor;
5. MEDLL, EML, and Double-Delta observation-level time.

The end-to-end real-time factor is
`processed_signal_seconds / wall_clock_seconds`.  The retained Python
benchmark in `scripts/baselines/runtime_benchmark.py` is a reference for
algorithmic scaling only; it is not a GPU result.  GPU numbers must come from
the actual receiver log or a timed CUDA run.

## 10. Result status and non-overclaim rule

At the time this protocol was added, the new baselines and the 560-cell
generalization grid had **no completed numerical result** in the public
repository.  The existing conditional PVT results remain valid, but they must
not be rewritten as if EML, Double-Delta, or MEDLL had already been measured.
