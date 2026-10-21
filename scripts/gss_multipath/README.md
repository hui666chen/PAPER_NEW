# Four-round GS multipath estimator — offline reproduction

Source: Xin Qi and Bing Xu, IEEE TIM 75 (2026), 8504218,
DOI 10.1109/TIM.2026.3666007, Section IV, equations 26–31 and Algorithms 1–2.
This is an independently written experimental implementation, not author code.

## Run

From the PAPER_NEW repository root:

```bash
python3 -m unittest scripts.gss_multipath.test_gss -v
python3 -m scripts.gss_multipath.run_scans --input /absolute/epoch_responses.csv --output /absolute/new_gss_results.json
```

The adapter reuses `read_epoch_scans` from the existing baseline module.
It accepts window-averaged files too, but those are not a reproduction of the
paper's coherent-epoch tracking-loop experiment. No simulator truth is consumed.
The default CSV axis is the receiver's local-code shift. It is reversed into
physical delay because the CUDA kernel uses code(t+shift). Use
`--axis physical-delay` only for files already using code(t-delay).

## Implemented algorithm

- Choose I/Q by absolute prompt amplitude; orient by the selected prompt sign.
- Maximize the selected signed correlation over [-1, 1.5] chips.
- Four searches with tolerance 0.001 chip; evaluate the current objective at
  peak+1 to resolve destructive components (26,29).
- Progressively subtract reference components with coefficients 1/4, 2/4, 3/4.
- Reject consecutive candidates separated by approximately one chip
  (strict 0.995–1.005), replacing the latest candidate by its predecessor,
  rolling back and using the squared coefficient in equation 31.
- Return the minimum candidate coordinate; preserve all intermediate results.

## Explicit transcription decisions

Algorithm 2 initializes f^1 but its first search line names f^j with j=2.
We search f^(j-1), then construct f^j, consistently with equations 27–29.
Algorithm 1 omits retained-point coordinate assignments while copying values;
we update both coordinates and values as required by golden-section search.
On a rejected candidate, its predecessor's amplitude is retained for a possible
subsequent rollback. The paper does not explicitly specify that stored amplitude.
These interpretations must be compared with author code if available.

## Limits and receiver integration

Default reference is a unit triangular L1 C/A ACF. `estimate` accepts a
normalized real reference callable for a later receiver-response ablation;
that change must be identified separately from the original reproduction.

The callable correlation interface permits future on-demand GPU correlations.
The CSV adapter uses linear interpolation and REJECTS insufficient measured
support: the paper's search plus peak+1 requires [-1,2.5] chips.
It never zero-pads missing measured samples or silently shrinks the search.
Interpolation cannot create sub-sample information; 401 taps alone do not
guarantee adequate support or resolution.

`build_scanner.py` builds a separate copy with receiver shifts [-3,+1] chips
(physical delays [-1,+3]); it does not modify the production binary.
`replay.py` and `score_pvt.py` perform the retained-data exploratory comparison.
See PVT_REPLAY_DESIGN.md for the fixed calibration and scoring conventions.

Output is a peak coordinate in the input delay axis. Its sign relative to the
receiver NCO/pseudorange, zero reference and units need independent calibration
before applying corrections. Numerical convergence is NOT physical confidence.
No existing gate, receiver configuration, tracking loop or PVT is modified.

Tests cover deterministic, ideal-response numerical behavior only. They do
not establish noisy multipath accuracy, real-time performance or PVT benefit.
The paper's multi-round heuristic is not guaranteed to find the true direct
path for arbitrary multipath or NLOS; failures must be retained in evaluation.
