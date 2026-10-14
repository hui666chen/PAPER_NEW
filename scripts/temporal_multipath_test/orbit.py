"""Minimal GPS broadcast-ephemeris tools (RINEX 3 navigation, IS-GPS-200 orbit).

Only what the temporal multipath test needs: satellite elevation/azimuth seen
from a static receiver.  Light-time and relativistic terms are ignored; the
resulting elevation error (well below 0.01 deg) is irrelevant here.
"""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path

import numpy as np

MU = 3.986005e14
OMEGA_E = 7.2921151467e-5
GPS_EPOCH = datetime(1980, 1, 6, tzinfo=timezone.utc)
LEAP_SECONDS_2017 = 18  # valid for dates after 2017-01-01


def utc_to_gps(t: datetime) -> tuple[int, float]:
    if t.tzinfo is None:
        t = t.replace(tzinfo=timezone.utc)
    sec = (t - GPS_EPOCH).total_seconds() + LEAP_SECONDS_2017
    week = int(sec // 604800)
    return week, sec - week * 604800


def _f(text: str) -> float:
    text = text.strip().replace("D", "E").replace("d", "e")
    return float(text) if text else 0.0


def read_gps_nav(path: Path) -> dict[int, list[dict]]:
    """Return {prn: [ephemeris dict, ...]} from a RINEX 3 (mixed or GPS) nav file."""
    lines = Path(path).read_text(errors="ignore").splitlines()
    start = next(i for i, ln in enumerate(lines) if "END OF HEADER" in ln) + 1
    eph: dict[int, list[dict]] = {}
    i = start
    while i < len(lines):
        ln = lines[i]
        if not ln.strip():
            i += 1
            continue
        system = ln[0]
        n_extra = {"G": 7, "E": 7, "C": 7, "J": 7, "I": 7, "R": 3, "S": 3}.get(system, 7)
        if system != "G":
            i += 1 + n_extra
            continue
        block = lines[i:i + 8]
        i += 8
        prn = int(block[0][1:3])
        vals = [_f(block[0][23 + 19 * k: 42 + 19 * k]) for k in range(3)]
        for row in block[1:]:
            vals += [_f(row[4 + 19 * k: 23 + 19 * k]) for k in range(4)]
        names = ["af0", "af1", "af2", "iode", "crs", "dn", "m0", "cuc", "e", "cus", "sqrta",
                 "toe", "cic", "omega0", "cis", "i0", "crc", "omega", "omegadot", "idot",
                 "l2codes", "week", "l2p", "acc", "health", "tgd", "iodc", "ttm", "fit"]
        rec = dict(zip(names, vals))
        if rec.get("health", 0) == 0 and rec["sqrta"] > 0:
            eph.setdefault(prn, []).append(rec)
    return eph


def sat_ecef(rec: dict, week: int, tow: np.ndarray) -> np.ndarray:
    t = np.asarray(tow, dtype=float) + (week - rec["week"]) * 604800.0
    tk = t - rec["toe"]
    a = rec["sqrta"] ** 2
    n = np.sqrt(MU / a ** 3) + rec["dn"]
    m = rec["m0"] + n * tk
    e = rec["e"]
    ek = m.copy()
    for _ in range(12):
        ek = m + e * np.sin(ek)
    nu = np.arctan2(np.sqrt(1 - e * e) * np.sin(ek), np.cos(ek) - e)
    phi = nu + rec["omega"]
    du = rec["cus"] * np.sin(2 * phi) + rec["cuc"] * np.cos(2 * phi)
    dr = rec["crs"] * np.sin(2 * phi) + rec["crc"] * np.cos(2 * phi)
    di = rec["cis"] * np.sin(2 * phi) + rec["cic"] * np.cos(2 * phi)
    u = phi + du
    r = a * (1 - e * np.cos(ek)) + dr
    inc = rec["i0"] + di + rec["idot"] * tk
    xp, yp = r * np.cos(u), r * np.sin(u)
    om = rec["omega0"] + (rec["omegadot"] - OMEGA_E) * tk - OMEGA_E * rec["toe"]
    return np.column_stack((xp * np.cos(om) - yp * np.cos(inc) * np.sin(om),
                            xp * np.sin(om) + yp * np.cos(inc) * np.cos(om),
                            yp * np.sin(inc)))


def lla_to_ecef(lat_deg: float, lon_deg: float, h: float) -> np.ndarray:
    a, f = 6378137.0, 1 / 298.257223563
    e2 = f * (2 - f)
    lat, lon = np.deg2rad(lat_deg), np.deg2rad(lon_deg)
    n = a / np.sqrt(1 - e2 * np.sin(lat) ** 2)
    return np.array([(n + h) * np.cos(lat) * np.cos(lon), (n + h) * np.cos(lat) * np.sin(lon),
                     (n * (1 - e2) + h) * np.sin(lat)])


def az_el(sat: np.ndarray, rx: np.ndarray, lat_deg: float, lon_deg: float) -> tuple[np.ndarray, np.ndarray]:
    lat, lon = np.deg2rad(lat_deg), np.deg2rad(lon_deg)
    d = sat - rx[None, :]
    e = -np.sin(lon) * d[:, 0] + np.cos(lon) * d[:, 1]
    n = (-np.sin(lat) * np.cos(lon) * d[:, 0] - np.sin(lat) * np.sin(lon) * d[:, 1]
         + np.cos(lat) * d[:, 2])
    u = np.cos(lat) * np.cos(lon) * d[:, 0] + np.cos(lat) * np.sin(lon) * d[:, 1] + np.sin(lat) * d[:, 2]
    el = np.rad2deg(np.arctan2(u, np.hypot(e, n)))
    az = np.rad2deg(np.arctan2(e, n)) % 360.0
    return az, el


def trajectories(nav: Path, start_utc: datetime, duration_s: int, lat: float, lon: float, alt: float,
                 min_el: float) -> dict[int, tuple[np.ndarray, np.ndarray]]:
    """{prn: (az_deg[t], el_deg[t])} for satellites above ``min_el`` for the whole span."""
    eph = read_gps_nav(nav)
    week, tow0 = utc_to_gps(start_utc)
    tow = tow0 + np.arange(duration_s, dtype=float)
    rx = lla_to_ecef(lat, lon, alt)
    out = {}
    for prn, recs in sorted(eph.items()):
        mid = tow0 + duration_s / 2
        rec = min(recs, key=lambda r: abs((week - r["week"]) * 604800 + mid - r["toe"]))
        az, el = az_el(sat_ecef(rec, week, tow), rx, lat, lon)
        if el.min() >= min_el:
            out[prn] = (az, el)
    return out


def synthetic_trajectories(duration_s: int, min_el: float, count: int = 4, seed: int = 1):
    """Fallback passes with a realistic 0.1-0.5 deg/min elevation rate (self-test only)."""
    rng = np.random.default_rng(seed)
    out = {}
    t = np.arange(duration_s) / 60.0
    for k in range(count):
        el0 = rng.uniform(min_el + 2, 70)
        rate = rng.uniform(0.1, 0.5) * rng.choice([-1, 1])
        el = np.clip(el0 + rate * t, min_el, 88)
        out[100 + k] = (np.full(duration_s, rng.uniform(0, 360)), el)
    return out
