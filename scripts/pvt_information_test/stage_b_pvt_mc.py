#!/usr/bin/env python3
"""Stage B: multi-satellite PVT Monte Carlo on Stage A cases (no IQ).

Question: given the per-satellite cases from Stage A, does putting the
profile-derived DD error interval into the PVT ("interval-constrained WLS")
beat the strongest classical alternatives that use the SAME curve
information, when several satellites carry short-delay multipath?

Per trial: one geometry (N satellites), K of them get an injected Stage A
case, the rest get clean Stage A cases.  Every satellite's error is
    e_i = e_DD_i (from Stage A) + n_i,   n_i ~ N(0, sigma_extra^2 (1 + 1/sin^2 el) / 2),
and every method sees the same e_i.  All satellites carry an interval,
including the clean ones (the method does not know which are affected).

Methods
  DD            WLS on DD pseudoranges
  EML           WLS on EML pseudoranges (e_EML from Stage A)
  DD_FDE        DD + chi-square RAIM fault exclusion (Pfa 1e-3)
  DD_DETEX      DD, exclude satellites whose curve detector fires (T > threshold)
  DD_DETW       DD, flagged satellites de-weighted (sigma x 5)
  MIDPOINT      DD corrected by the interval midpoint
  INTERVAL      proposed: WLS with a bounded bias per satellite, b_i in [lo_i, hi_i]
  INTERVAL_DET  secondary: bounded bias only on flagged satellites
  ORACLE        true error removed (lower bound, not a method)

Frozen decision (see protocols/SHORT_DELAY_PVT_INFORMATION_TEST_PROTOCOL.md):
  primary scenario = N 9, K 3, random placement, delays {0.05, 0.10} chip,
  amplitudes {0.25, 0.50}, sigma_extra 1.0 m.
  PASS if INTERVAL 3D RMSE <= 0.85 x best(DD, EML, DD_FDE, DD_DETEX, DD_DETW)
  with the upper 95% bootstrap bound of that ratio < 1.0, AND in the K=0
  scenario INTERVAL <= 1.05 x DD.  Otherwise STOP the interval-PVT direction.

Usage
    python3 scripts/pvt_information_test/stage_b_pvt_mc.py \
        --cases results/pvt_information_test/stage_a_cases.csv
    optional real sky: --rtklib-stat path/to/solution.pos.stat   (rnx2rtkp -y 2)
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import lsq_linear
from scipy.stats import chi2

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
SEED = 2026101421
METHODS = ("DD", "EML", "DD_FDE", "DD_DETEX", "DD_DETW", "MIDPOINT", "INTERVAL",
           "INTERVAL_DET", "ORACLE")
COMPARATORS = ("DD", "EML", "DD_FDE", "DD_DETEX", "DD_DETW")
PRIMARY = {"n_sat": 9, "k": 3, "placement": "random", "delays": (0.05, 0.10),
           "amps": (0.25, 0.50), "sigma_extra": 1.0}
CASE_DTYPE = np.dtype([
    ("scene", "U16"), ("delay_chip", "f8"), ("amplitude", "f8"),
    ("e_dd_m", "f8"), ("e_eml_m", "f8"),
    ("dd_interval_lo_m", "f8"), ("dd_interval_hi_m", "f8"),
    ("flagged", "?"), ("T", "f8"),
])


def read_cases(path: Path) -> np.ndarray:
    """Read Stage A rows without requiring pandas in the experiment environment."""
    rows = []
    with Path(path).open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            rows.append((row["scene"], float(row["delay_chip"]), float(row["amplitude"]),
                         float(row["e_dd_m"]), float(row["e_eml_m"]),
                         float(row["dd_interval_lo_m"]), float(row["dd_interval_hi_m"]),
                         row["flagged"].strip().lower() == "true", float(row["T"])))
    if not rows:
        raise ValueError(f"no Stage A rows in {path}")
    return np.asarray(rows, dtype=CASE_DTYPE)


# ------------------------------------------------------------------ geometry
def random_geometries(count: int, n_sat: int, rng, max_pdop: float = 3.5, min_el: float = 10.0):
    out = []
    while len(out) < count:
        az = rng.uniform(0, 360, n_sat)
        s = rng.uniform(np.sin(np.deg2rad(min_el)), 1.0, n_sat)
        el = np.rad2deg(np.arcsin(s))
        h = design(az, el)
        try:
            pdop = np.sqrt(np.trace(np.linalg.inv(h.T @ h)[:3, :3]))
        except np.linalg.LinAlgError:
            continue
        if pdop <= max_pdop:
            out.append((az, el))
    return out


def geometries_from_rtklib_stat(path: Path, min_el: float = 10.0, every: int = 10):
    """Read $SAT lines of an RTKLIB solution-status file (rnx2rtkp -y 2)."""
    epochs: dict[float, list[tuple[float, float]]] = {}
    for line in Path(path).read_text(errors="ignore").splitlines():
        if not line.startswith("$SAT"):
            continue
        p = line.split(",")
        # $SAT,week,tow,sat,frq,az,el,resp,resc,vsat,snr,...
        tow, az, el, vsat = float(p[2]), float(p[5]), float(p[6]), int(p[9])
        if vsat and el >= min_el:
            epochs.setdefault(tow, []).append((az, el))
    keys = sorted(epochs)[::every]
    geoms = [(np.array([a for a, _ in epochs[k]]), np.array([e for _, e in epochs[k]]))
             for k in keys if len(epochs[k]) >= 6]
    if not geoms:
        raise ValueError(f"no usable $SAT epochs in {path}")
    return geoms


def design(az_deg, el_deg):
    az, el = np.deg2rad(az_deg), np.deg2rad(el_deg)
    los = np.column_stack((np.cos(el) * np.sin(az), np.cos(el) * np.cos(az), np.sin(el)))
    return np.column_stack((-los, np.ones(len(az))))


# ------------------------------------------------------------------ solvers
def wls(h, z, w):
    hw = h * w[:, None]
    return np.linalg.solve(h.T @ hw, hw.T @ z)


def wls_fde(h, z, w, pfa=1e-3):
    keep = np.ones(len(z), bool)
    while True:
        hk, zk, wk = h[keep], z[keep], w[keep]
        x = wls(hk, zk, wk)
        r = zk - hk @ x
        dof = keep.sum() - 4
        if dof < 1 or float(np.sum(wk * r * r)) <= chi2.ppf(1 - pfa, dof) or keep.sum() <= 5:
            return x
        s = np.linalg.inv(hk.T @ (hk * wk[:, None]))
        lev = np.einsum("ij,jk,ik->i", hk, s, hk) * wk
        nr = np.abs(r) * np.sqrt(wk) / np.sqrt(np.maximum(1 - lev, 1e-9))
        idx = np.flatnonzero(keep)[int(np.argmax(nr))]
        keep[idx] = False


def wls_exclude(h, z, w, flagged, score, min_keep=5):
    keep = ~flagged
    if keep.sum() < min_keep:
        order = np.argsort(score)  # keep the least suspicious satellites
        keep = np.zeros(len(z), bool)
        keep[order[:min_keep]] = True
    return wls(h[keep], z[keep], w[keep])


def wls_interval(h, z, w, lo, hi, active=None, eps=1e-3):
    """min sum w (z - Hx - b)^2 + eps sum w b^2,  lo <= b <= hi on active satellites."""
    n = len(z)
    act = np.ones(n, bool) if active is None else active
    m = int(act.sum())
    if m == 0:
        return wls(h, z, w)
    sw = np.sqrt(w)
    s = np.zeros((n, m))
    s[np.flatnonzero(act), np.arange(m)] = 1.0
    top = np.hstack((h * sw[:, None], s * sw[:, None]))
    reg = np.hstack((np.zeros((m, 4)), np.diag(np.sqrt(eps) * sw[act])))
    a = np.vstack((top, reg))
    rhs = np.concatenate((z * sw, np.zeros(m)))
    lb = np.concatenate((np.full(4, -np.inf), lo[act]))
    ub = np.concatenate((np.full(4, np.inf), hi[act]))
    ub = np.maximum(ub, lb + 1e-9)
    sol = lsq_linear(a, rhs, bounds=(lb, ub), method="bvls")
    return sol.x[:4]


# ------------------------------------------------------------------ Monte Carlo
def run_scenario(cases: np.ndarray, geoms, sc: dict, trials: int, rng) -> dict:
    clean = cases[cases["scene"] == "clean"]
    inj = cases[(cases["scene"] == "injected") & np.isin(cases["delay_chip"], sc["delays"])
                & np.isin(cases["amplitude"], sc["amps"])]
    err3 = {m: np.empty(trials) for m in METHODS}
    errh = {m: np.empty(trials) for m in METHODS}
    errv = {m: np.empty(trials) for m in METHODS}
    for t in range(trials):
        az, el = geoms[rng.integers(len(geoms))]
        n = min(sc["n_sat"], len(az))
        pick = rng.choice(len(az), n, replace=False)
        az, el = az[pick], el[pick]
        h = design(az, el)
        # RTKLIB-style elevation model: sigma^2 = s^2 (1 + 1/sin^2 el) / 2 (= s at zenith)
        sig = sc["sigma_extra"] * np.sqrt(0.5 * (1.0 + 1.0 / np.sin(np.deg2rad(el)) ** 2))
        w = 1.0 / sig ** 2
        k = min(sc["k"], n - 5)
        if sc["placement"] == "lowest":
            aff = np.argsort(el)[:k]
        else:
            aff = rng.choice(n, k, replace=False)
        rows = clean[rng.integers(len(clean), size=n)].copy()
        if k > 0:
            rows[aff] = inj[rng.integers(len(inj), size=k)]
        noise = rng.normal(0.0, sig)
        e_dd = rows["e_dd_m"]
        e_eml = rows["e_eml_m"]
        lo = rows["dd_interval_lo_m"]
        hi = rows["dd_interval_hi_m"]
        flag = rows["flagged"]
        score = rows["T"]
        z = e_dd + noise
        sols = {
            "DD": wls(h, z, w),
            "EML": wls(h, e_eml + noise, w),
            "DD_FDE": wls_fde(h, z, w),
            "DD_DETEX": wls_exclude(h, z, w, flag, score),
            "DD_DETW": wls(h, z, np.where(flag, w / 25.0, w)),
            "MIDPOINT": wls(h, z - 0.5 * (lo + hi), w),
            "INTERVAL": wls_interval(h, z, w, lo, hi),
            "INTERVAL_DET": wls_interval(h, z, w, lo, hi, active=flag),
            "ORACLE": wls(h, noise, w),
        }
        for m, x in sols.items():
            err3[m][t] = np.linalg.norm(x[:3])
            errh[m][t] = np.linalg.norm(x[:2])
            errv[m][t] = abs(x[2])
    res = {"scenario": sc, "trials": trials, "methods": {}}
    for m in METHODS:
        res["methods"][m] = {
            "rmse_3d_m": float(np.sqrt(np.mean(err3[m] ** 2))),
            "rmse_h_m": float(np.sqrt(np.mean(errh[m] ** 2))),
            "rmse_v_m": float(np.sqrt(np.mean(errv[m] ** 2))),
            "p95_3d_m": float(np.quantile(err3[m], 0.95)),
        }
    best = min(COMPARATORS, key=lambda m: res["methods"][m]["rmse_3d_m"])
    boot = []
    for _ in range(2000):
        i = rng.integers(trials, size=trials)
        boot.append(np.sqrt(np.mean(err3["INTERVAL"][i] ** 2)) / np.sqrt(np.mean(err3[best][i] ** 2)))
    res["best_comparator"] = best
    res["interval_over_best_ratio"] = res["methods"]["INTERVAL"]["rmse_3d_m"] / res["methods"][best]["rmse_3d_m"]
    res["interval_over_best_ratio_ci95"] = [float(np.quantile(boot, 0.025)), float(np.quantile(boot, 0.975))]
    res["interval_over_dd_ratio"] = res["methods"]["INTERVAL"]["rmse_3d_m"] / res["methods"]["DD"]["rmse_3d_m"]
    return res


def scenarios():
    yield dict(PRIMARY, name="primary_K3")
    yield dict(PRIMARY, k=0, name="clean_K0")
    for k in (1, 2, 4):
        yield dict(PRIMARY, k=k, name=f"K{k}")
    yield dict(PRIMARY, placement="lowest", name="primary_K3_lowest_elevation")
    for d in (0.02, 0.05, 0.10, 0.20):
        yield dict(PRIMARY, delays=(d,), amps=(0.10, 0.25, 0.50), name=f"K3_delay_{d:.2f}")
    for s in (0.5, 2.0):
        yield dict(PRIMARY, sigma_extra=s, name=f"primary_K3_sigma_{s}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cases", type=Path, default=REPO / "results/pvt_information_test/stage_a_cases.csv")
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--trials", type=int, default=4000)
    ap.add_argument("--geometries", type=int, default=200)
    ap.add_argument("--rtklib-stat", type=Path, default=None)
    ap.add_argument("--allow-stage-a-stop", action="store_true",
                    help="run even if Stage A said STOP (result is then exploratory only)")
    args = ap.parse_args()
    out = args.out or args.cases.parent
    summary_path = args.cases.parent / "stage_a_summary.json"
    stage_a = json.loads(summary_path.read_text())["summary"]["stage_a_decision"] if summary_path.exists() else "unknown"
    if stage_a != "STAGE_A_PASS" and not args.allow_stage_a_stop:
        sys.exit(f"Stage A decision is {stage_a}; frozen rule says stop. "
                 f"Use --allow-stage-a-stop only for an explicitly exploratory run.")
    rng = np.random.default_rng(SEED)
    cases = read_cases(args.cases)
    geoms = (geometries_from_rtklib_stat(args.rtklib_stat) if args.rtklib_stat
             else random_geometries(args.geometries, PRIMARY["n_sat"], rng))
    results = []
    for sc in scenarios():
        r = run_scenario(cases, geoms, sc, args.trials, rng)
        results.append(r)
        mm = r["methods"]
        print(f"[stage B] {sc['name']:<30} " + " ".join(
            f"{m}={mm[m]['rmse_3d_m']:.2f}" for m in METHODS), flush=True)
    prim = next(r for r in results if r["scenario"]["name"] == "primary_K3")
    k0 = next(r for r in results if r["scenario"]["name"] == "clean_K0")
    gain_ok = prim["interval_over_best_ratio"] <= 0.85 and prim["interval_over_best_ratio_ci95"][1] < 1.0
    safe_ok = k0["interval_over_dd_ratio"] <= 1.05
    decision = "STAGE_B_PASS" if (gain_ok and safe_ok) else "STAGE_B_STOP_INTERVAL_PVT"
    payload = {"status": "complete", "stage_a_decision": stage_a, "decision": decision,
               "rule": "primary INTERVAL/best_comparator <= 0.85 with CI upper < 1.0, and K0 INTERVAL/DD <= 1.05",
               "geometry_source": str(args.rtklib_stat) if args.rtklib_stat else "random_pdop_le_3.5",
               "exploratory_only": stage_a != "STAGE_A_PASS", "results": results}
    out.mkdir(parents=True, exist_ok=True)
    (out / "stage_b_pvt_summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    lines = ["# Stage B: multi-satellite PVT Monte Carlo", "", f"Decision: **{decision}**  ",
             f"Rule: {payload['rule']}  ", f"Stage A: {stage_a}; geometry: {payload['geometry_source']}", "",
             "3D RMSE (m)", "", "| scenario | " + " | ".join(METHODS) + " | best comparator | INTERVAL/best [95% CI] |",
             "|---|" + "---:|" * (len(METHODS) + 2)]
    for r in results:
        mm = r["methods"]
        ci = r["interval_over_best_ratio_ci95"]
        lines.append(f"| {r['scenario']['name']} | " + " | ".join(f"{mm[m]['rmse_3d_m']:.2f}" for m in METHODS)
                     + f" | {r['best_comparator']} | {r['interval_over_best_ratio']:.3f} [{ci[0]:.3f}, {ci[1]:.3f}] |")
    lines += ["", "Semi-synthetic: per-satellite errors and intervals come from Stage A "
              "(PRN09 clean windows plus known injections); geometry and extra noise are simulated. "
              "This decides whether a hardware-in-the-loop experiment is worth running; it is not itself a PVT result."]
    (out / "STAGE_B_REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"[stage B] decision: {decision} -> {out}")


if __name__ == "__main__":
    main()
