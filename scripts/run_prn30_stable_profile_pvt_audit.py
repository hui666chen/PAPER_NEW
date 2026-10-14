#!/usr/bin/env python3
"""Post-hoc PRN30 PVT sensitivity to stable four-column profiling.

Uses the already retained common single-path RINEX as the immutable baseline,
changes only G30 C1C using the B-stable profile estimates, and writes every
new file under PAPER_NEW/results. Existing receiver/PVT artifacts are read-only.
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys


PAPER = Path(__file__).resolve().parents[2]
REPO = PAPER / "PAPER_NEW"
ROOT = PAPER / "gps_l1_prn30_formal_pvt"
COMMON = ROOT / "final_sign_persistence/common_observable_pvt"
OUT = REPO / "results/shared_ratio_numerical_audit_20261014/pvt_stable_B_corrected"
sys.path.insert(0, str(ROOT))

import build_common_observable_rinex as rinex  # noqa: E402
import score_common_observable_pvt as pvt_score  # noqa: E402


SCENES = {
    "clean": range(60, 76),
    "amp_0p01": range(79, 95),
    "amp_0p05": range(100, 116),
    "amp_0p10": range(121, 137),
}
RESCORE = REPO / "results/shared_ratio_numerical_audit_20261014/prn30_retained_four_column_rescore_corrected.csv"
SOURCE_RINEX = COMMON / "observations/common_single.26O"
NAV = COMMON / "navigation/GSDR269k11.26N"
CONFIG = COMMON / "rnx2rtkp_single.conf"
SOLVER = COMMON / "bin/rnx2rtkp"
CHIP_METRES = 299_792_458.0 / 1.023e6


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_stable_schedules() -> dict[str, dict[int, float]]:
    rows = list(csv.DictReader(RESCORE.open(newline="", encoding="utf-8")))
    stable = {
        (r["scene"], int(r["window"])): r
        for r in rows if r["solver"] == "B_stable_with_affine"
    }
    schedules = {"cross_only_stable_B": {}, "sign_persistent_stable_B": {}}
    for scene, seconds in SCENES.items():
        history = []
        for window in seconds:
            row = stable[(scene, window)]
            correction = float(row["correction_m"])
            cross_active = row["cross_active"] == "True"
            schedules["cross_only_stable_B"][window] = correction if cross_active else 0.0
            history.append((window, row, correction, cross_active))
            history = history[-4:]
            sign = 0 if abs(correction) < 1e-15 else (1 if correction > 0 else -1)
            persistent = False
            if len(history) == 4 and cross_active and sign:
                agreeing = sum(
                    active and ((1 if value > 0 else -1) if abs(value) >= 1e-15 else 0) == sign
                    for _, _, value, active in history
                )
                persistent = agreeing >= 3
            schedules["sign_persistent_stable_B"][window] = (
                correction if persistent else 0.0
            )
    return schedules


def create_variant(source: Path, destination: Path,
                   schedule: dict[int, float]) -> dict:
    lines = source.read_text(encoding="ascii").splitlines(keepends=True)
    output = []
    in_header = True
    second = None
    changed = 0
    for line in lines:
        if in_header:
            output.append(line)
            if "END OF HEADER" in line:
                in_header = False
            continue
        if line.startswith(">"):
            second = rinex.relative_second(rinex.parse_epoch(line))
            output.append(line)
            continue
        if line.startswith("G30"):
            if second is None:
                raise ValueError("G30 observation before first epoch")
            original = float(line[3:17])
            correction = schedule.get(second, 0.0)
            adjusted = original + correction
            replaced = line[:3] + rinex.format_c1c(adjusted, line[17:19]) + line[19:]
            if replaced[19:] != line[19:]:
                raise AssertionError("non-C1C RINEX fields changed")
            if abs(correction) > 0:
                changed += 1
            output.append(replaced)
        else:
            output.append(line)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text("".join(output), encoding="ascii")
    return {"target_rows_changed": changed, "output_sha256": sha256(destination)}


def run_solver(observation: Path, solution: Path, log: Path) -> None:
    solution.parent.mkdir(parents=True, exist_ok=True)
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("w", encoding="utf-8") as stream:
        subprocess.run(
            [str(SOLVER), "-k", str(CONFIG), "-o", str(solution),
             str(observation), str(NAV)],
            cwd=ROOT, check=True, stdout=stream, stderr=subprocess.STDOUT,
        )


def scene_metric(solution_path: Path, scene: str) -> dict:
    rows = pvt_score.parse_solution(solution_path)
    selected = [rows[second] for second in SCENES[scene] if second in rows]
    if len(selected) != 16:
        raise RuntimeError(f"{solution_path}: {scene} has {len(selected)}/16 epochs")
    return pvt_score.metric_row(scene, solution_path.stem, selected)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    schedules = load_stable_schedules()
    source_hash = sha256(SOURCE_RINEX)
    artifacts = {}
    for method, schedule in schedules.items():
        obs = OUT / "observations" / f"{method}.26O"
        solution = OUT / "solutions" / f"{method}.pos"
        log = OUT / "logs" / f"{method}.log"
        artifacts[method] = create_variant(SOURCE_RINEX, obs, schedule)
        run_solver(obs, solution, log)
        artifacts[method].update({
            "observation_path": str(obs),
            "solution_path": str(solution),
            "solution_sha256": sha256(solution),
            "log_path": str(log),
        })

    original_solutions = {
        "single_A_baseline": COMMON / "solutions/single.pos",
        "cross_only_A": COMMON / "solutions/cross_only_refined.pos",
        "sign_persistent_A": COMMON / "solutions/sign_persistent_refined.pos",
        "cross_only_stable_B": OUT / "solutions/cross_only_stable_B.pos",
        "sign_persistent_stable_B": OUT / "solutions/sign_persistent_stable_B.pos",
    }
    metrics = []
    for scene in SCENES:
        for method, path in original_solutions.items():
            if "stable_B" in method:
                continue
            metrics.append(scene_metric(path, scene) | {"method": method})
        for method in schedules:
            metrics.append(scene_metric(OUT / "solutions" / f"{method}.pos", scene)
                           | {"method": method})
    metric_path = OUT / "pvt_metrics.csv"
    with metric_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(metrics[0]))
        writer.writeheader()
        writer.writerows(metrics)

    by_scene = {}
    for scene in SCENES:
        rows = {r["method"]: r for r in metrics if r["scene"] == scene}
        baseline = float(rows["single_A_baseline"]["rmse_3d_m"])
        by_scene[scene] = {
            method: {
                "rmse_3d_m": float(row["rmse_3d_m"]),
                "improvement_vs_single_m": baseline - float(row["rmse_3d_m"]),
            }
            for method, row in rows.items()
        }
    manifest = {
        "schema": 1,
        "status": "posthoc_PRN30_PVT_sensitivity_audit",
        "confirmatory": False,
        "changed_only": "G30 C1C in copies of the retained common single-path RINEX",
        "source_rinex": str(SOURCE_RINEX),
        "source_rinex_sha256": source_hash,
        "response_rescore_csv": str(RESCORE),
        "response_rescore_sha256": sha256(RESCORE),
        "config_sha256": sha256(CONFIG),
        "navigation_sha256": sha256(NAV),
        "solver_sha256": sha256(SOLVER),
        "scenes_seconds": {k: [v.start, v.stop - 1] for k, v in SCENES.items()},
        "fixed_gate": "P_cross <= 0.20; causal same-sign 3-of-4 persistence for persistent variant",
        "variants": artifacts,
        "metrics": by_scene,
        "caution": (
            "This is a post-hoc sensitivity comparison on one retained PRN30 "
            "noise realization, not an independent confirmation. No original "
            "files were modified and no PRN9 claim is adjudicated."
        ),
    }
    manifest_path = OUT / "audit_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2), flush=True)


if __name__ == "__main__":
    main()
