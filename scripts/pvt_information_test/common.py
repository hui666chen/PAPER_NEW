"""Shared pieces for the short-delay PVT information test.

Delay convention (same as the existing PAPER_NEW scripts):
  * the 401 taps are an open-loop scan centred on the receiver prompt;
  * a direct path whose tap position is ``tau`` is the response ``model.exact(-tau)``;
  * a reflection with separation ``delta`` is ``model.exact(delta - tau)``;
  * the Double-Delta / EML ``tau_chip`` lives on the same tap axis, so for a
    DD-corrected pseudorange the residual error is
        e_DD = (tau_direct - tau_DD + c_DD) * CHIP_METRES,
    where ``c_DD`` is the clean-calibration median of (tau_DD - tau_onepath).

The two-path profile below is a self-contained copy of the stable
variable-projection profile in
``scripts/short_delay_detectability/refine_bias_intervals.py`` (orthogonal QR
projection of the constant/linear tap nuisance, no raw-Gram pseudoinverse).
It is copied, not imported, so this test does not depend on the external
``shared_ratio_existence`` module.
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

CHIP_METRES = 299_792_458.0 / 1_023_000.0
TAP_COUNT = 401
SEPARATIONS = np.round(np.arange(0.04, 1.2000001, 0.02), 10)
RHO_MAGNITUDES = (0.01, 0.02, 0.05, 0.10, 0.20, 0.34, 0.50)
RHO_PHASES = np.arange(16) * (2.0 * np.pi / 16.0)
RHO_GRID = np.asarray(
    [0.0] + [m * np.exp(1j * p) for m in RHO_MAGNITUDES for p in RHO_PHASES],
    dtype=np.complex128,
)


def quantile_higher(values, q: float) -> float:
    """NumPy-version-compatible ``higher`` quantile (same rule on both APIs)."""
    try:
        return float(np.quantile(values, q, method="higher"))
    except TypeError:  # NumPy < 1.22 used ``interpolation`` for this option.
        return float(np.quantile(values, q, interpolation="higher"))


def affine_basis(tap: np.ndarray) -> np.ndarray:
    one = np.ones(len(tap), dtype=np.complex128)
    q, _ = np.linalg.qr(np.column_stack((one, tap.astype(np.complex128))), mode="reduced")
    return q


def project_out(y: np.ndarray, q: np.ndarray) -> np.ndarray:
    # Column-vector projection: y - Q (Q^H y).  Valid for complex Q.
    return y - q @ (q.conj().T @ y)


class StableProfile:
    """Two-path RSS profiled over complex gain, shared ratio and separation.

    ``profile(y)`` returns, on the direct-delay grid ``self.tau``:
      * ``two``: min RSS over (separation, rho) at each tau (rho=0 included);
      * ``one``: single-path RSS at each tau.
    """

    def __init__(self, model, tap: np.ndarray, step: float = 0.002,
                 tau_lo: float = -0.40, tau_hi: float = 0.40,
                 separations: np.ndarray = SEPARATIONS,
                 rho_grid: np.ndarray = RHO_GRID):
        self.model = model
        self.step = float(step)
        self.tau = np.round(np.arange(tau_lo, tau_hi + step * 0.5, step), 12)
        self.separation = np.asarray(separations, dtype=float)
        self.rho = np.asarray(rho_grid, dtype=np.complex128)
        self.q = affine_basis(tap)
        args = np.concatenate((
            -self.tau,
            (self.separation[None, :] - self.tau[:, None]).ravel(),
        ))
        self.arguments = np.unique(np.round(args, 10))
        index = {round(float(v), 10): i for i, v in enumerate(self.arguments)}
        response = np.asarray(model.exact(self.arguments)).T  # (taps, n_args)
        projected = response - self.q @ (self.q.conj().T @ response)
        self.templates = projected
        self.gram = projected.conj().T @ projected
        self.norm = np.maximum(np.diag(self.gram).real, 1e-30)
        self.d = np.asarray([index[round(float(-t), 10)] for t in self.tau])
        self.r = np.asarray([[index[round(float(s - t), 10)] for s in self.separation]
                             for t in self.tau])
        self.cross = self.gram[self.d[:, None], self.r]

    def profile(self, y: np.ndarray) -> dict[str, np.ndarray | float]:
        y_perp = project_out(np.asarray(y, dtype=np.complex128), self.q)
        corr = self.templates.conj().T @ y_perp
        energy = float(np.vdot(y_perp, y_perp).real)
        gd = self.norm[self.d]
        gr = self.norm[self.r]
        vd = corr[self.d]
        vr = corr[self.r]
        one = np.maximum(energy - np.abs(vd) ** 2 / gd, 0.0)
        two = one.copy()
        chunk = 16
        for start in range(0, len(self.rho), chunk):
            rc = self.rho[start:start + chunk][:, None, None]
            den = (gd[None, :, None] + rc * self.cross[None] + np.conj(rc) * self.cross.conj()[None]
                   + (np.abs(rc) ** 2) * gr[None])
            num = np.abs(vd[None, :, None] + np.conj(rc) * vr[None]) ** 2
            cost = np.maximum(energy - num / np.maximum(den.real, 1e-30), 0.0)
            two = np.minimum(two, cost.min(axis=(0, 2)))
        return {"one": one, "two": two, "energy": energy, "corr": corr}

    def onepath_fit(self, prof: dict) -> tuple[float, complex]:
        """Single-path delay (parabolic refinement on the grid) and complex gain."""
        one = prof["one"]
        k = int(np.argmin(one))
        tau = float(self.tau[k])
        if 0 < k < len(one) - 1:
            a, b, c = one[k - 1], one[k], one[k + 1]
            den = a - 2.0 * b + c
            if den > 0:
                tau += 0.5 * self.step * (a - c) / den
        gain = complex(prof["corr"][self.d[k]] / self.norm[self.d[k]])
        return tau, gain

    def interval(self, prof: dict, kappa: float) -> tuple[float, float, float]:
        two = prof["two"]
        m = float(two.min())
        keep = self.tau[two <= m * (1.0 + kappa) + 1e-14 * max(m, 1.0)]
        return float(keep.min()), float(keep.max()), float(self.tau[int(np.argmin(two))])

    def required_kappa(self, prof: dict, target_tau: float) -> float:
        """Smallest kappa whose retained-tau hull contains ``target_tau``."""
        two = prof["two"]
        m = max(float(two.min()), 1e-30)
        half = 0.5 * self.step
        left = two[self.tau <= target_tau + half]
        right = two[self.tau >= target_tau - half]
        if not left.size or not right.size:
            return float("inf")
        return max(0.0, max(float(left.min()), float(right.min())) / m - 1.0)

    @staticmethod
    def detection_stat(prof: dict) -> float:
        return float(prof["one"].min() / max(float(prof["two"].min()), 1e-30) - 1.0)


def read_epoch_matrix(path: Path, wanted_windows, tap_count: int = TAP_COUNT) -> dict[int, np.ndarray]:
    """Same CSV reader as refine_bias_intervals.py (columns window, epoch, tap, mixed_i, mixed_q)."""
    wanted = set(int(w) for w in wanted_windows)
    rows: dict[int, dict[int, dict[int, complex]]] = {w: {} for w in wanted}
    with Path(path).open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            window = int(row["window"])
            if window not in wanted:
                continue
            rows[window].setdefault(int(row["epoch"]), {})[int(row["tap"])] = complex(
                float(row["mixed_i"]), float(row["mixed_q"]))
    out = {}
    for window in sorted(wanted):
        epochs = sorted(rows[window])
        if not epochs:
            raise ValueError(f"{path}: window {window} missing")
        mat = []
        for e in epochs:
            taps = rows[window][e]
            if len(taps) != tap_count:
                raise ValueError(f"{path}: window {window} epoch {e} has {len(taps)} taps")
            mat.append([taps[t] for t in range(tap_count)])
        out[window] = np.asarray(mat, dtype=np.complex128)
    return out


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
