# Huber + binary gate: exploratory retained-data experiment

No IQ generation; reuse matching original 401-tap, 100-epoch responses from
seed 2026101922. Existing results and receiver binary remain unchanged.

Model: original shared complex reflection ratio, per-epoch complex gain.
Use explicitly gain-only stable projection, matching the old effective model
whose unscaled pseudoinverse discarded affine columns. Check LS agreement
against saved fits before making algorithm-comparison claims.

Search unchanged: tau -0.40:0.02:0.40 chip (H0 uses 0.01), separation
0.04:0.02:1.20 chip; original 113 complex rho candidates, up to 0.50.
Do NOT truncate the search at 0.40: fit first, reject afterwards.

Complex radial Huber: L(r)=|r|² for |r|<=delta; otherwise
2 delta |r|-delta². Delta=1.345 times component MAD scale of training H0
LS residual, with floor 1e-12. Shared delta for H0/H1 and held-out scoring.
Full grid joint IRLS: at most 6 steps, relative objective convergence 1e-5.
This is an iterative approximate optimizer, not certified global Huber.

Even/odd cross fitting unchanged. Inherited gain threshold
0.0032589351980997073 is exploratory under Huber, not a calibrated false-alarm
claim. Single-window evidence requires both held-out gains above threshold,
no delay boundary, and converged fits. Rho >0.40 in either fitted fold is a
hard veto, both on evidence and after 3-of-5 persistence. Otherwise correction
is zero. No scene-label/oracle amplitude veto and no three-level decisions.

Report LS/Huber ordinary RSS and robust losses separately; losses of different
types must not be compared numerically as if they were identical. Report delay
boundary, rho veto, convergence, evidence and activation counts. A true 0.50
scene may still activate if rho is underestimated; this is a failure to be
reported, not suppressed.

Calibrate clean zero with windows 27..39 only; evaluate PVT on 40..79 with
identical observations and settings. Clean must actually traverse the gate.
Existing data have been inspected: all outcomes are exploratory, not a new
independent confirmation. Do not claim PVT benefit before computing it.

Paper direction is provisional robust fitting + binary safety gating.
GSS remains a failed exploratory implementation / future investigation,
not a demonstrated positive method.
