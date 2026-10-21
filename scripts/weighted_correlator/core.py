"""Coherent weighted early/late primitive, not a full MAA-DLL reproduction.

Implements the weighted discriminator structure of Bengherabi et al. (2025),
IEEE Access, DOI 10.1109/ACCESS.2025.3581856, eq. 12. Weights and full E-L
spacings are explicit inputs; no BOC lookup table is applied to BPSK C/A.
Input axis MUST be physical delay in chips, increasing to the right.
"""
from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class Configuration:
    spacings: tuple = (0.1, 0.2)
    weights: tuple = (1.0, -0.5)
    radius: float = 0.15
    step: float = 0.001

    def __post_init__(self):
        if (len(self.spacings) != len(self.weights) or not self.spacings
                or not np.all(np.isfinite(self.spacings + self.weights))
                or min(self.spacings) <= 0
                or not np.any(np.asarray(self.weights) != 0)
                or not np.isfinite(self.radius + self.step)
                or not 0 < self.step < self.radius):
            raise ValueError("invalid spacings, weights, or search support")


def estimate(axis, correlation, *, reference=0.0, config=Configuration()):
    """Find nearest positive-slope discriminator root to tracking reference.

Carrier wipe-off uses the observed reference prompt, NOT injected phase.
Missing support / flat / absent roots return invalid with no correction.
The result is a local code coordinate, not a calibrated pseudorange bias.
No truth, clean data, path count or path parameter is used here.
"""
    x, z = np.asarray(axis, float), np.asarray(correlation, complex)
    if (x.ndim != 1 or x.shape != z.shape or len(x) < 5
            or not np.all(np.isfinite(x)) or not np.all(np.isfinite(z))
            or np.any(np.diff(x) <= 0) or not np.isfinite(reference)):
        raise ValueError("invalid correlation samples")
    result = dict(valid=False, estimate_chip=None, root_count=0, reason="")
    margin = max(config.spacings) / 2
    if (reference-config.radius-margin < x[0]
            or reference+config.radius+margin > x[-1]):
        return dict(result, reason="insufficient_measured_support")
    scale = np.max(np.abs(z))
    if scale == 0:
        return dict(result, reason="zero_signal")
    z = z/scale
    prompt = np.interp(reference, x, z)
    if abs(prompt) < 1e-10:
        return dict(result, reason="undefined_prompt_phase")
    coherent = np.real(z*np.conj(prompt)/abs(prompt))
    grid = np.linspace(reference-config.radius, reference+config.radius,
                       int(np.ceil(2*config.radius/config.step))+1)
    values = sum(w*(np.interp(grid-d/2, x, coherent)
                    - np.interp(grid+d/2, x, coherent))
                 for d, w in zip(config.spacings, config.weights))
    # Reject zero plateaus: they carry no timing information. Positive slope
    # corresponds to E-L on an isolated positive correlation peak.
    slope = np.diff(values)/np.diff(grid)
    roots = []
    for i in range(len(grid)-1):
        if values[i] <= 0 <= values[i+1] and slope[i] > 1e-8:
            root = grid[i] - values[i]/slope[i]
            if grid[0] < root < grid[-1] and (not roots or abs(root-roots[-1]) > 1e-9):
                roots.append(float(root))
    if not roots:
        return dict(result, reason="no_stable_local_root")
    return dict(valid=True, estimate_chip=min(roots, key=lambda r: abs(r-reference)),
                root_count=len(roots), reason="local_root_not_confidence")
