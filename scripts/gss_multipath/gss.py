"""Offline Qi & Xu (TIM 2026, DOI 10.1109/TIM.2026.3666007) reproduction.

Equations 26--31; four adaptive subtraction rounds. Code phase is expressed
in the supplied correlation-delay coordinate, NOT a signed pseudorange fix.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Callable
import math

import numpy as np


ScalarFunction = Callable[[float], float]


def triangular_acf(delay: float) -> float:
    """Unit-height ideal GPS L1 C/A reference; receiver filtering is not modeled."""
    return max(0.0, 1.0 - abs(delay))


@dataclass(frozen=True)
class Config:
    lower: float = -1.0
    upper: float = 1.5
    tolerance: float = 0.001
    rejection_tolerance: float = 0.005
    max_iterations: int = 128

    def validate(self):
        values = (self.lower, self.upper, self.tolerance, self.rejection_tolerance)
        if not all(math.isfinite(v) for v in values):
            raise ValueError("nonfinite configuration")
        if not self.lower < self.upper or self.tolerance <= 0:
            raise ValueError("invalid search interval or tolerance")
        if self.rejection_tolerance <= 0 or self.max_iterations < 1:
            raise ValueError("invalid rejection tolerance or iteration limit")


@dataclass
class Peak:
    delay: float
    value: float
    lower: float
    upper: float
    evaluations: int
    converged: bool


def golden_maximum(function: ScalarFunction, config: Config) -> Peak:
    """Algorithm 1 with retained point coordinates updated alongside values.

    Convergence certifies interval shrinkage only, not a global/physical peak.
    """
    config.validate()
    count = 0

    def evaluate(x):
        nonlocal count
        count += 1
        result = float(function(float(x)))
        if not math.isfinite(result):
            raise ValueError("nonfinite objective")
        return result

    a, b = config.lower, config.upper
    inverse_phi = (math.sqrt(5.0) - 1.0) / 2.0
    x1, x2 = b - inverse_phi * (b - a), a + inverse_phi * (b - a)
    f1, f2 = evaluate(x1), evaluate(x2)
    for _ in range(config.max_iterations):
        if b - a <= config.tolerance:
            break
        if f1 >= f2:
            b, x2, f2 = x2, x1, f1
            x1 = b - inverse_phi * (b - a)
            f1 = evaluate(x1)
        else:
            a, x1, f1 = x1, x2, f2
            x2 = a + inverse_phi * (b - a)
            f2 = evaluate(x2)
    delay = (a + b) / 2
    return Peak(delay, evaluate(delay), a, b, count, b - a <= config.tolerance)


def estimate(
    correlation: Callable[[float], complex],
    reference: ScalarFunction = triangular_acf,
    config: Config = Config(),
) -> dict:
    """Run four rounds on a callable (including future direct GPU correlation).

    reference must be real, centered at zero and normalized to unit height.
    No truth, multipath amplitude, clean counterfactual or old gate is read.
    """
    config.validate()
    if not np.isclose(reference(0.0), 1.0, atol=1e-8):
        raise ValueError("reference must have unit value at zero")
    prompt = complex(correlation(0.0))
    if not np.isfinite(prompt) or abs(prompt) == 0:
        raise ValueError("zero/nonfinite prompt cannot select an arm and sign")
    arm = "I" if abs(prompt.real) >= abs(prompt.imag) else "Q"
    sign = float(np.sign(prompt.real if arm == "I" else prompt.imag))

    def initial(t):
        value = complex(correlation(t))
        return sign * (value.real if arm == "I" else value.imag)

    objectives = [initial]
    candidates, amplitudes, trace = [], [], []
    for r in range(4):
        objective = objectives[r]
        peak = golden_maximum(objective, config)
        if not peak.converged:
            raise ValueError("golden search exceeded iteration limit")
        raw_delay = peak.delay
        alias = bool(r and abs(abs(raw_delay - candidates[r - 1]) - 1.0)
                     < config.rejection_tolerance)
        shifted = False
        next_value = None
        if alias:
            # Algorithm 2 lines 13--18: replace the invalid candidate by the
            # preceding estimate. Eq. 31 rolls back to the prior objective.
            candidate, amplitude = candidates[r - 1], amplitudes[r - 1]
            if r < 3:
                base = objectives[r - 1]
                coefficient = ((r + 1) / 4.0) ** 2 * amplitude
        else:
            # Eq. 26/29, evaluated on the current (not original) objective.
            next_value = abs(objective(raw_delay + 1.0))
            shifted = peak.value < next_value
            candidate = raw_delay + 1.0 if shifted else raw_delay
            amplitude = -next_value if shifted else peak.value
            if r < 3:
                base = objective
                coefficient = (r + 1) / 4.0 * amplitude
        candidates.append(float(candidate))
        amplitudes.append(float(amplitude))
        trace.append({
            "round": r + 1, "search": asdict(peak),
            "candidate_chip": float(candidate), "amplitude": float(amplitude),
            "plus_one_absolute_value": next_value,
            "shifted_plus_one": bool(shifted), "one_chip_alias_rejected": alias,
        })
        if r < 3:
            def updated(t, base=base, coefficient=coefficient, candidate=candidate):
                return base(t) - coefficient * reference(t - candidate)
            objectives.append(updated)
    return {
        "estimate_chip": min(candidates), "arm": arm, "sign": sign,
        "rounds": trace, "config": asdict(config),
        "numerically_converged": True,
        "physical_confidence_calibrated": False,
        "pseudorange_correction_applied": False,
    }


class SampledCorrelation:
    """Linear interpolation with strict support; never silently extrapolates."""
    def __init__(self, delay_chip, values):
        self.x = np.asarray(delay_chip, dtype=float)
        self.y = np.asarray(values, dtype=complex)
        if (self.x.ndim != 1 or self.x.size < 5 or self.y.shape != self.x.shape
                or not np.all(np.isfinite(self.x)) or not np.all(np.isfinite(self.y))
                or np.any(np.diff(self.x) <= 0)):
            raise ValueError("invalid sampled correlation")

    def __call__(self, delay):
        if delay < self.x[0] - 1e-12 or delay > self.x[-1] + 1e-12:
            raise ValueError(f"query {delay:.6f} outside measured support "
                             f"[{self.x[0]}, {self.x[-1]}] chip")
        return complex(np.interp(delay, self.x, self.y))

    def validate_support(self, config):
        # Search [-1,1.5] plus the +1 query needs measured support to +2.5.
        required_min = min(0.0, config.lower)
        required_max = max(0.0, config.upper + 1.0)
        if self.x[0] > required_min + 1e-12 or self.x[-1] < required_max - 1e-12:
            raise ValueError(f"paper search requires [{required_min}, {required_max}] "
                             f"chip; scan supplies [{self.x[0]}, {self.x[-1]}]")
