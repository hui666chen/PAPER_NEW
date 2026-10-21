# Initial implementation check (2026-10-20)

This is an implementation diagnostic, not an experiment demonstrating PVT gain.

Command: `python3 -m unittest scripts.gss_multipath.test_gss -v`

Coverage: golden maximum against analytic maxima; clean shifts, signs and
I/Q arm selection; short-delay phase sweep; three weak reflections; destructive
peak relocation; one-chip over-subtraction rejection; sampled support checking;
CSV adapter; explicit nonfinite/zero/nonconvergent failures.

Counterexample retained: ideal triangular direct path at 0 plus a reflection
at +0.1 chip with amplitude 1.2 and phase 180 degrees returns approximately
+0.100084 chip, selecting the reflection. The prompt sign is reversed in this
case. This is not a passing mitigation result and needs comparison with the
authors' implementation/assumptions before making reproduction claims.

Existing input inspected:
`results/prn09_phase_amp_delay_grid_20261019/seed_2026101922/pvt_gated_shared_ratio_d0p05_a0p20_p090_hostgpu/clean/correlator_scan/chunk_00/frozen_open_loop_prn9_responses.csv`.
It is a 100-epoch window-average format on [-2,+2] chips. The strict adapter
rejects this support for the original [-1,+1.5] search with +1 queries.
Need wider correlation outputs or on-demand correlation from retained IQ.
No new IQ generation is required merely to compute those extra correlations.

No changes were made to the existing gate, receiver binary, IQ, PVT outputs,
or historical experiment reports.
