# Weighted early/late discriminator: initial port

Status: tested numerical primitive only; NOT deployed, NOT a validated PVT
improvement, and NOT a reproduction of the full 2025 MAA-DLL algorithm.

Reference: Bengherabi et al., Efficient Multipath Mitigation in GNSS Using
Dynamic Selection and Weighting of Early-Late Correlators, IEEE Access 2025,
https://doi.org/10.1109/ACCESS.2025.3581856, equation 12 (EDD structure).
The EDD structure predates this paper. No novelty is claimed for this port.

The configurable sum of coherent early-minus-late terms accepts measured
complex correlation samples. Full E-L spacings and signed coefficients are
explicit. Default (0.1,0.2) chips with (1,-0.5) coefficients is classical DD,
not an optimized three-pair EDD or the paper's adaptive MAA-DLL.
It must be compared separately from the repository's power-domain DD with
0.1/0.5 chip spacings. Neither may silently replace historical baselines.

The paper's full method first estimates an echo delay via residual correlation,
then selects weights/spacing/subcarrier through a BOC lookup table. That table
cannot be transplanted unchanged into GPS L1 C/A. Adaptive selection and
three-pair weight optimization are NOT implemented here. No coefficients
have been selected against retained test outcomes.

Run: `python3 -m unittest scripts.weighted_correlator.test_core -v`

Input axis is physical delay; reverse and negate CUDA code(t+shift) scan axes.
Outputs are code coordinates, NOT correction metres. Independent clean-zero
calibration and sign verification are required before applying to pseudorange.
Invalid/unsupported/flat cases produce no estimate. A found root is not a
confidence guarantee, especially under destructive interference or false lock.
The local +/-0.15 chip support is a prototype setting, not tuned on PVT.
Interpolating taps does not add measurement information. Front-end bandwidth,
DLL dynamics, data-bit handling and carrier lock remain integration issues.

## Practical receiver option already present in the GPU receiver

The GPS L1 GPU tracking block already has a noncoherent E-L DLL discriminator
and configurable tap shifts. Its `early_late_space_chips` setting is the
one-sided tap offset: the actual E and L taps are at `-value` and `+value`.
Therefore a 0.10-chip total E-L spacing corresponds to setting
`Tracking_1C.early_late_space_chips=0.05`. The checked-in receiver config
currently uses 0.5, i.e. a 1.0-chip total spacing. A one-line config patch is
provided as `narrow_eml_receiver.patch`; it is not applied to the receiver tree.

In the retained one-seed offline PVT comparison, the 0.10-chip EML baseline
reduced 3D RMSE by 2.9% to 17.7% versus the uncorrected receiver across the
four injected-reflection cases, while the fixed weighted candidate did not
beat EML in any of the five scenes. This supports using narrow EML as the first
practical receiver setting to validate. It is not a closed-loop deployment
result: final acceptance still requires tracking/PVT on an independent run and
checking clean performance. Original Huber/model gates were not changed.
