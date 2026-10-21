"""Post-hoc round-by-round clean diagnostic; never alters frozen predictions."""
import csv
import itertools
import json

import numpy as np
from scripts.gss_multipath.gss import SampledCorrelation, estimate
from scripts.gss_multipath.replay import OUT, METRES


def main():
    records, examples = [], []
    for source in sorted((OUT/"clean").glob("chunk_*/epochs.csv")):
        with source.open() as stream:
            for key, group in itertools.groupby(csv.DictReader(stream),
                    lambda r: (int(r["window"]), int(r["epoch"]))):
                rows = list(group)
                x = -np.array([float(r["delay_chip"]) for r in rows])[::-1]
                z = np.array([complex(float(r["mixed_i"]),float(r["mixed_q"]))
                              for r in rows])[::-1]
                result = estimate(SampledCorrelation(x, z))
                peaks = [r["candidate_chip"] for r in result["rounds"]]
                records.append([*key, *peaks, result["estimate_chip"]])
                if len(examples) < 3 and peaks[0]-result["estimate_chip"] > .3:
                    examples.append({"window":key[0],"epoch":key[1],**result})
    a = np.array(records)
    summary = {
        "epochs": len(a),
        "first_peak_mean_chip": float(a[:,2].mean()),
        "final_mean_chip":float(a[:,-1].mean()),
        "first_peak_sd_chip":float(a[:,2].std()),
        "final_sd_chip":float(a[:,-1].std()),
        "later_round_moves_earlier_gt_0p05_chip":int(np.sum(a[:,2]-a[:,-1]>.05)),
        "later_round_moves_earlier_gt_0p1_chip":int(np.sum(a[:,2]-a[:,-1]>.1)),
        "interpretation": "Later subtraction rounds introduce early candidates in clean data. "
            "This localizes the failure; it does not distinguish interpolation, reference shape, "
            "noise or paper transcription effects.",
        "examples": examples,
    }
    with (OUT/"CLEAN_ROUND_DIAGNOSTIC.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["window","epoch","round1","round2","round3","round4","selected"])
        writer.writerows(records)
    (OUT/"CLEAN_ROUND_DIAGNOSTIC.json").write_text(json.dumps(summary,indent=2))
    print(json.dumps({k:v for k,v in summary.items() if k!="examples"},indent=2))


if __name__=="__main__":
    main()
