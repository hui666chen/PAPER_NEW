# Exploratory replay design, 2026-10-20

Fixed before reading GSS replay performance: seed 2026101922, clean and four
retained PRN9 scenes (delay .05, amplitude .2 with phases 0/90/180, and amplitude
.5 phase90). Windows27–39 calibrate each method's constant clean zero by mean.
Windows40–79 evaluate, subject to source RINEX availability.

Same GPU scanner, same 401-tap samples, same 100 coherent epochs per window
for GS, existing EML and DD. Arithmetic mean of epoch estimates per window.
No fitting of phase/amplitude, no new threshold or safety gate. All methods
change only G09 C1C; all other observations, ephemerides and RTKLIB settings
remain identical within each scene. Original receiver C1C is the fourth branch.

The scanner's local code is code(t+shift), whereas the paper's axis is code(t-delay).
Reverse/reorder the scanner axis. Existing baseline_schedule uses -shift*c/fcode,
so the corresponding physical-delay correction is +delay*c/fcode.
Subtract the clean calibrated zero before this mapping. Do not select correction
sign based on PVT performance.

Report availability and metrics on the common solved seconds across all methods.
No independent-run confidence interval is possible with one seed. These data
were used for previous diagnostics; this is not a fresh confirmatory experiment.
Any paired clean-curve displacement is a counterfactual proxy, not exact
geometric pseudorange truth. Ground-truth position is used only after predictions
and PVT solver outputs are fixed.

This test uses interpolation on an offline replay with the existing scanner's
front-end approximation, not on-demand correlations inside the DLL. It cannot
establish full reproduction of the paper's tracking-loop or real-time results.
