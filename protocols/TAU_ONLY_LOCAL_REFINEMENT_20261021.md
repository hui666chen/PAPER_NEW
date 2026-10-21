# Tau-only local refinement: retained-data diagnostic

Primary estimator is the existing stable gain-only shared-ratio LS model.
No Huber/GSS method changes, no new IQ and no PVT run in this diagnostic.
Use all five retained scenes, windows27..79; clean27..39 calibrates the median
zero, clean40..79 evaluates zero scatter. Compare all scenes' fold disagreement
on40..79. No validation labels enter estimation or tuning.

Each fold starts at its saved LS grid winner. Hold separation and the full
complex rho fixed. Reprofile only the existing per-epoch complex gain for
each tau. Evaluate the SAME bank's `exact`, not a different interpolation model.
H1 local interval is tau_grid +/-0.01 chip, clipped to original [-0.40,0.40].
For a fair H0/H1 test, similarly refine H0 tau within +/-0.005 chip, its half-grid
cell. Original candidate and both endpoints remain candidates; bounded scalar
minimization uses xatol1e-7 chip and maxiter100. Never accept higher training RSS.
No global search-domain expansion, no separation/rho/phase reoptimization.

Gate code, CV threshold,3-of-5,estimated-rho cap0.40 and correction sign are
unchanged. Recompute held-out RSS/CV and global boundary flags from refined
parameters; local interval-edge counts are diagnostics, not a new veto.
Report recalculated gate decisions but do not tune them or run PVT in this step.

Report clean std before/after, median and RMS two-fold disagreement, training
RSS improvement, held-out changes, local edge rate and gate activation. No
invented pass threshold: reductions are descriptive on one inspected seed.
Conditional tau improvement is not evidence that joint identifiability is fixed.
Original results remain immutable; new outputs use a separate directory.
