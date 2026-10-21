# Retained-IQ Huber / binary gate diagnostic

Run from repository root:

```
OPENBLAS_NUM_THREADS=1 python3 -m unittest discover -s scripts/huber_binary_gate -p test_fit.py -v
OPENBLAS_NUM_THREADS=1 python3 scripts/huber_binary_gate/run.py
python3 scripts/huber_binary_gate/score.py
```

This reuses the original receiver's retained per-epoch correlator responses
from the existing IQ runs. It does not generate IQ, change receiver tracking,
or deploy Huber into the real-time GPU receiver. The offline robust fitting
uses CPU BLAS; no GPU execution claim is made.

`run.py` can be restricted with `--scenes` and `--windows`. Existing per-window
outputs are skipped; for changed estimator settings use a new `--out`
directory. Do not mix different parameter versions in one directory.
`score.py` requires all five scenarios and all 53 windows.

No affine nuisance columns are silently reintroduced. Stable gain-only LS
is compared against the old saved effective pseudoinverse fits in the audit.
The robust cost is complex radial Huber, not independent clipping of I and Q.
Complete parameter-grid search is repeated during joint IRLS. It is still
possible to converge to a local Huber solution; no global-optimality claim.

Output includes an uncapped Huber ablation. Only `huber_cap` is the requested
binary candidate. Original persistence is retained, and the fitted-amplitude
veto is applied after persistence as well. Parameter cap is not ground truth.

The inherited gain threshold is not independently calibrated for this new
loss. These reused-data diagnostics cannot establish deployment safety,
publication novelty, significance, or superiority over DD.
