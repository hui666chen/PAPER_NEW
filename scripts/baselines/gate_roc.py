#!/usr/bin/env python3
"""Evaluate single-window and K-of-M gate operating points.

The score is observation-only.  Scene labels are accepted only for the
post-hoc evaluation table and are never used to create a decision column.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def _truth(row: dict[str, str], positive_scenes: set[str]) -> bool:
    return row.get("scene", "") in positive_scenes


def _decision(row: dict[str, str], threshold: float) -> bool:
    score = float(row["cv_gain_min"])
    no_boundary = row.get("cv_no_boundary", "true").lower() in {"1", "true", "yes"}
    return bool(score > threshold and no_boundary)


def _rates(rows: list[dict[str, object]]) -> dict[str, float | int]:
    positives = [row for row in rows if row["truth"]]
    negatives = [row for row in rows if not row["truth"]]
    tp = sum(bool(row["decision"]) for row in positives)
    fp = sum(bool(row["decision"]) for row in negatives)
    return {
        "tp": tp,
        "fp": fp,
        "positive_count": len(positives),
        "negative_count": len(negatives),
        "tpr": tp / len(positives) if positives else float("nan"),
        "fpr": fp / len(negatives) if negatives else float("nan"),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--records-csv", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--positive-scenes", default="amp_0p05,amp_0p10")
    parser.add_argument("--score-field", default="cv_gain_min")
    args = parser.parse_args()

    positive_scenes = {value for value in args.positive_scenes.split(",") if value}
    with args.records_csv.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError("empty records CSV")
    if args.score_field != "cv_gain_min":
        raise ValueError("the frozen gate score is cv_gain_min")
    scores = sorted({float(row[args.score_field]) for row in rows})
    thresholds = [scores[0] - 1.0] + [
        (left + right) / 2.0 for left, right in zip(scores, scores[1:])
    ] + [scores[-1] + 1.0]

    single_window = []
    for threshold in thresholds:
        decisions = [
            {"truth": _truth(row, positive_scenes), "decision": _decision(row, threshold)}
            for row in rows
        ]
        single_window.append({"threshold": threshold, **_rates(decisions)})

    # K-of-M is evaluated within each seed/scene sequence.  If a seed column
    # is absent (old compact records), one sequence per scene is used.
    sequences: dict[tuple[str, str], list[dict[str, str]]] = {}
    for row in rows:
        key = (row.get("seed", "single_sequence"), row.get("scene", "unknown"))
        sequences.setdefault(key, []).append(row)
    for sequence in sequences.values():
        sequence.sort(key=lambda row: int(float(row.get("window", 0))))

    persistence = []
    for n in (1, 3, 5, 7):
        for k in range(1, n + 1):
            for threshold in thresholds:
                decisions = []
                for sequence in sequences.values():
                    raw = [_decision(row, threshold) for row in sequence]
                    for index, row in enumerate(sequence):
                        active = index >= n - 1 and sum(raw[index - n + 1:index + 1]) >= k
                        decisions.append({
                            "truth": _truth(row, positive_scenes),
                            "decision": active,
                        })
                persistence.append({
                    "n": n,
                    "k": k,
                    "threshold": threshold,
                    **_rates(decisions),
                })

    result = {
        "truth_accessed_for_decision": False,
        "positive_scenes_for_scoring": sorted(positive_scenes),
        "score_field": args.score_field,
        "single_window": single_window,
        "k_of_m": persistence,
        "sqm_status": "not_available_unless_an_independent_sqm_column_is_supplied",
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "single_window_points": len(single_window),
        "k_of_m_points": len(persistence),
        "output": str(args.output_json),
    }, indent=2))


if __name__ == "__main__":
    main()
