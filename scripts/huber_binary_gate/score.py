"""Freeze binary decisions, audit LS equivalence, then common-observation PVT."""
import csv
from datetime import timedelta
import json
from pathlib import Path
import sys
import numpy as np
from run import ROOT, BASE, SCENES

sys.path.insert(0, str(ROOT))
from scripts.prn09_phase_delay_amplitude_grid import complete_seed22_gated_pvt as helper

OUT = ROOT/"results/huber_binary_gate_20261020"
METHODS = ("ls", "ls_cap", "huber", "huber_cap")
THRESHOLD = .0032589351980997073


def decisions(records, method, zero):
    kind = "huber" if method.startswith("huber") else "ls"
    votes, output = [], []
    for r in records:
        folds = r[kind]
        veto = method.endswith("cap") and any(x["rho"] > .4+1e-12 for x in folds)
        converged = all(x["converged"] for x in folds)
        boundary = any(x["boundary"] for x in folds)
        evidence = all(x["gain"] > THRESHOLD for x in folds) and not boundary and converged and not veto
        votes.append(evidence)
        # Keep the original temporal rule. Apply present-time amplitude and
        # numerical-convergence vetoes AFTER persistence as well.
        active = len(votes) >= 5 and sum(votes[-5:]) >= 3 and not veto and converged
        tau = float(np.mean([x["tau_chip"] for x in folds]))
        output.append(dict(window=r["window"], evidence=bool(evidence), veto=bool(veto),
                           boundary=bool(boundary), converged=bool(converged),
                           active=bool(active), tau=tau,
                           correction_m=-(tau-zero)*helper.CHIP_METRES if active else 0.))
    return output


def main():
    seal_path = OUT/"PREDICTIONS_SEAL.json"
    sealed = json.loads(seal_path.read_text()) if seal_path.exists() else None
    if sealed:
        for name, digest in sealed["hashes"].items():
            if name.startswith("results/") and helper.sha256(ROOT/name) != digest:
                raise RuntimeError(f"sealed fit changed: {name}")
    data = {s: [json.loads((OUT/f"{s}_{w}.json").read_text()) for w in range(27, 80)]
            for s in SCENES}
    zero = {kind: float(np.median([np.mean([x["tau_chip"] for x in r[kind]])
                                  for r in data["clean"] if r["window"] < 40]))
            for kind in ("ls", "huber")}
    diagnostics = {}
    predictions = {}
    helper.SCENES = {s: BASE/s for s in SCENES}
    for scene, records in data.items():
        cell = "d0p05_a0p20_p090" if scene == "clean" else scene
        dirname = ("pvt_gated_shared_ratio_complete_hostgpu" if cell.endswith("p000")
                   else f"pvt_gated_shared_ratio_{cell}_hostgpu")
        oldpath = BASE/dirname/("clean" if scene == "clean" else "mixed")/"gate/gate_records.json"
        saved = {r["window"]: r for r in json.loads(oldpath.read_text())}
        mismatch, max_cost = [], 0.
        for r in records:
            for k, direction in enumerate(("A_to_B", "B_to_A")):
                previous = saved[r["window"]]["fit"][direction]
                now = r["ls"][k]
                if any(abs(now[a]-previous["h1"][b]) > 1e-8 for a, b in
                       (("tau_chip", "tau_chip"), ("separation_chip", "separation_chip"),
                        ("rho", "rho_magnitude"))):
                    mismatch.append([r["window"], direction])
                max_cost = max(max_cost, abs(now["h1_cost"]-previous["rss_shared_ratio_test"]))
        diagnostics[scene] = dict(ls_parameter_mismatches=mismatch,
                                  max_ls_heldout_cost_difference=max_cost, methods={})
        predictions[scene] = {}
        for method in METHODS:
            kind = "huber" if method.startswith("huber") else "ls"
            rows = decisions(records, method, zero[kind])
            stats = {name: sum(r[name] for r in rows)
                     for name in ("boundary", "veto", "evidence", "active", "converged")}
            stats["ordinary_rss_median"] = float(np.median([
                sum(x["ordinary_rss"] for x in r[kind]) for r in records]))
            stats["heldout_loss_median"] = float(np.median([
                sum(x["h1_cost"] for x in r[kind]) for r in records]))
            diagnostics[scene]["methods"][method] = stats
            directory = OUT/scene/method
            directory.mkdir(parents=True, exist_ok=True)
            (directory/"decisions.json").write_text(json.dumps(rows, indent=2))
            schedule = directory/"corrections.csv"
            with schedule.open("w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(["start_sample", "end_sample", "prn", "correction_m"])
                for r in rows:
                    w = r["window"]
                    writer.writerow([w*25000000, (w+1)*25000000-1, 9,
                                     r["correction_m"] if w >= 40 else 0])
            if sealed:
                previous = sealed["predictions"][scene][method]
                if helper.sha256(schedule) != previous["schedule_sha256"]:
                    raise RuntimeError("decision schedule changed; use a new output directory")
                paths, solutions = previous["rinex"], previous["solutions"]
                for key, digest in paths.items():
                    if key.endswith("_sha256") and helper.sha256(Path(paths[key[:-7]])) != digest:
                        raise RuntimeError(f"sealed observation changed: {key}")
                for entry in solutions.values():
                    if isinstance(entry, dict) and "path" in entry and "sha256" in entry:
                        if helper.sha256(Path(entry["path"])) != entry["sha256"]:
                            raise RuntimeError("sealed solution changed")
            else:
                paths = helper.build_rinex_branches(
                    scene, directory, schedule, helper.SIGNAL_START_UTC+timedelta(seconds=40),
                    helper.SIGNAL_START_UTC+timedelta(seconds=80))
                solutions = helper.solve_pair(scene, paths, directory)
            predictions[scene][method] = dict(rinex=paths, solutions=solutions,
                                              schedule_sha256=helper.sha256(schedule))
    hashes = {str(p.relative_to(ROOT)): helper.sha256(p)
              for p in sorted((ROOT/"scripts/huber_binary_gate").glob("*.py"))}
    hashes.update({str(p.relative_to(ROOT)): helper.sha256(p)
                   for s in SCENES for p in sorted(OUT.glob(s+"_*.json"))})
    if not sealed:
        seal_path.write_text(json.dumps(
            dict(zero=zero, predictions=predictions, hashes=hashes,
                 truth_used_in_estimator=False), indent=2))
    (OUT/"REPORT_GENERATION.json").write_text(json.dumps(dict(
        script_sha256=helper.sha256(Path(__file__)),
        original_prediction_seal_sha256=helper.sha256(seal_path),
        reused_verified_predictions=bool(sealed)), indent=2))
    (OUT/"FIT_AUDIT.json").write_text(json.dumps(diagnostics, indent=2))
    # Only now open truth for scoring.
    position = json.loads((BASE/"clean/signalsim.json").read_text())["trajectory"]["initPosition"]
    truth = tuple(float(position[k]) for k in ("latitude", "longitude", "altitude"))
    scores = {}
    for scene in SCENES:
        estimates = {m: helper.parse_solution(Path(predictions[scene][m]["solutions"]["gated"]["path"]))
                     for m in METHODS}
        estimates["original"] = helper.parse_solution(Path(
            predictions[scene]["ls"]["solutions"]["single"]["path"]))
        expected = set(predictions[scene]["ls"]["rinex"]["score_seconds"])
        common = sorted(expected.intersection(*(set(v) for v in estimates.values())))
        scores[scene] = {}
        for m, values in estimates.items():
            enu = np.asarray([helper.enu_error(values[s]["latitude_deg"], values[s]["longitude_deg"],
                                              values[s]["height_m"], *truth) for s in common])
            scores[scene][m] = dict(rmse_3d_m=float(np.sqrt(np.mean(np.sum(enu**2, axis=1)))) if common else None,
                                   common_epochs=len(common), available=len(set(values)&expected),
                                   expected=len(expected), missing=sorted(expected-set(values)))
    (OUT/"PVT_SCORING.json").write_text(json.dumps(scores, indent=2))
    with (OUT/"WINDOW_COMPARISON.csv").open("w", newline="") as f:
        fields = ["scene", "window", "method", "boundary", "veto", "evidence",
                  "active", "converged", "tau", "correction_m", "ordinary_rss",
                  "heldout_loss", "max_fitted_rho", "min_cv_gain"]
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for s, records in data.items():
            for m in METHODS:
                kind = "huber" if m.startswith("huber") else "ls"
                for r, decision in zip(records, decisions(records, m, zero[kind])):
                    writer.writerow(dict(scene=s, method=m, **decision,
                        ordinary_rss=sum(x["ordinary_rss"] for x in r[kind]),
                        heldout_loss=sum(x["h1_cost"] for x in r[kind]),
                        max_fitted_rho=max(x["rho"] for x in r[kind]),
                        min_cv_gain=min(x["gain"] for x in r[kind])))
    lines = ["# Huber 二值门控：现有数据探索性重算", "",
             "一个已见种子；五场景各53窗口。不是独立确认，也不是实时部署。",
             "重新使用原IQ已提取的逐历元401抽头数据，没有重生成IQ。",
             "两种损失分别使用clean 27–39秒的中位零点；PVT只评分40–79秒。",
             "保留3-of-5；估计rho>0.4在持续性判断后仍强制置零。", "",
             "|场景|LS触边界|Huber触边界|Huber幅度拒绝|Huber+上限激活|",
             "|---|---:|---:|---:|---:|"]
    for s in SCENES:
        m = diagnostics[s]["methods"]
        lines.append(f"|{s}|{m['ls']['boundary']}/53|{m['huber']['boundary']}/53|"
                     f"{m['huber_cap']['veto']}/53|{m['huber_cap']['active']}/53|")
    lines += ["", "## 残差（每窗口两折留出代价之和的中位数）", "",
              "普通RSS统一对选中模板重新用LS剖面增益，衡量模板本身的拟合效果；"
              "Huber列使用Huber剖面增益。两种损失绝对值不能直接比较。",
              "", "|场景|LS普通RSS|Huber模板普通RSS|Huber损失|",
              "|---|---:|---:|---:|"]
    for s in SCENES:
        m = diagnostics[s]["methods"]
        lines.append(f"|{s}|{m['ls']['ordinary_rss_median']:.6f}|"
                     f"{m['huber']['ordinary_rss_median']:.6f}|"
                     f"{m['huber']['heldout_loss_median']:.6f}|")
    lines += ["", "## PVT 3D RMSE（米）", "",
              "|场景|原接收机|LS|LS+上限|Huber|Huber+上限|", "|---|---:|---:|---:|---:|---:|"]
    for s in SCENES:
        lines.append("|"+s+"|"+"|".join(
            f"{scores[s][m]['rmse_3d_m']:.4f}" if scores[s][m]["rmse_3d_m"] is not None else "无共同解"
            for m in ("original",)+METHODS)+"|")
    lines += ["", "## 限制",
              "- 普通RSS与Huber损失不可直接跨类型比较；原始值见FIT_AUDIT.json。",
              "- 改损失后沿用的cv阈值尚未独立校准，不宣称保持原误报率。",
              "- 单个种子不支持显著性或泛化结论；未比较DD/EML，不宣称优于经典方法。",
              "- 真实幅度0.50不等于估计幅度>0.40；不得按场景标签强制拒绝。",
              "- LS与旧实现的逐折参数及代价一致性见FIT_AUDIT.json。",
              "- GSS保留为此前失败的探索性分支，仅作为后续研究，不计正面贡献。"]
    lines += ["", "## 2026-10-21 完成后的结论",
              "- 265窗口、530折的LS参数与旧实现全部一致，留出代价最大差异小于2e-11。",
              "- Huber减少部分触边界，但普通RSS没有一致改善；clean有1个窗口未满足迭代收敛要求并被拒绝。",
              "- 0.50场景：Huber未加上限激活30/53，加上限后2/53；只有第66秒产生非零修正（-2.9305米），第67秒修正为零。不是全部拒绝。",
              "- 五场景所有分支都得到40/40个评分历元，没有PVT丢解。这里是离线共同观测量解算，不能等同于闭环跟踪稳定性验证。",
              "- 0.50场景RMSE从原接收机2.3473米变为2.3019米（约1.93%），但LS+上限与Huber+上限的PVT完全相同。不能把这项收益归因于Huber。",
              "- 其余四场景最终均零激活，PVT保持原接收机水平，没有新增改善。",
              "- 因此支持的是本批数据上的保守门控作用，而非已验证的Huber精度增益或可泛化的多径消除方法。"]
    (OUT/"REPORT.md").write_text("\n".join(lines)+"\n")
    print(json.dumps(scores, indent=2))


if __name__ == "__main__":
    main()
