"""Common-observation PVT diagnostic; no loop-feedback or improvement claims."""
import csv
from datetime import timedelta
import json
from pathlib import Path

import numpy as np
from scripts.gss_multipath.replay import OUT, SEED, SCENES, METRES
from scripts.prn09_phase_delay_amplitude_grid import complete_seed22_gated_pvt as helper


def main():
    clean = json.loads((OUT/"clean/window_estimates.json").read_text())
    methods = ("gss", "eml", "dd")
    zero = {method: float(np.mean([r[f"mixed_{method}"] for r in clean
                                  if r["window"] < 40])) for method in methods}
    calibration = {
        "zero_chip": zero, "calibration_windows": list(range(27, 40)),
        "heldout_windows": list(range(40, 80)),
        "mapping": "C1C += (physical_delay - clean_zero) * c/fcode",
        "sign_basis": "CUDA code(t+shift); physical_delay=-shift; existing baseline_schedule applies -shift",
        "scope": "same retained seed, disjoint clean time windows; exploratory",
    }
    (OUT/"PVT_CALIBRATION.json").write_text(json.dumps(calibration, indent=2))
    helper.SCENES = {scene: SEED/scene for scene in SCENES}
    predictions = {}
    for scene in SCENES:
        rows = json.loads((OUT/scene/"window_estimates.json").read_text())
        predictions[scene] = {}
        for method in methods:
            directory = OUT/scene/f"pvt_{method}"
            directory.mkdir(exist_ok=True)
            manifest_path = directory/"manifest.json"
            if manifest_path.exists():
                predictions[scene][method] = json.loads(manifest_path.read_text())
                continue
            schedule = directory/"corrections.csv"
            with schedule.open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=[
                    "start_sample", "end_sample", "prn", "correction_m"])
                writer.writeheader()
                for r in rows:
                    w = int(r["window"])
                    correction = (r[f"mixed_{method}"]-zero[method])*METRES if w >= 40 else 0
                    writer.writerow({"start_sample": w*25000000,
                                     "end_sample": (w+1)*25000000-1,
                                     "prn": 9, "correction_m": correction})
            paths = helper.build_rinex_branches(
                scene, directory, schedule,
                helper.SIGNAL_START_UTC+timedelta(seconds=40),
                helper.SIGNAL_START_UTC+timedelta(seconds=80))
            solutions = helper.solve_pair(scene, paths, directory)
            predictions[scene][method] = {"rinex": paths, "solutions": solutions,
                                         "schedule_sha256": helper.sha256(schedule)}
            manifest_path.write_text(json.dumps(predictions[scene][method], indent=2))
    (OUT/"PVT_PREDICTIONS_SEAL.json").write_text(json.dumps({
        "calibration": calibration, "predictions": predictions,
        "truth_used_in_estimator": False,
    }, indent=2))
    # Scoring only, after all outputs are fixed.
    init = json.loads((SEED/"clean/signalsim.json").read_text())["trajectory"]["initPosition"]
    truth = tuple(float(init[k]) for k in ("latitude", "longitude", "altitude"))
    report, errors = {}, []
    for scene in SCENES:
        estimates = {
            method: helper.parse_solution(Path(predictions[scene][method]["solutions"]["gated"]["path"]))
            for method in methods}
        estimates["original"] = helper.parse_solution(Path(
            predictions[scene]["gss"]["solutions"]["single"]["path"]))
        expected = set(predictions[scene]["gss"]["rinex"]["score_seconds"])
        common = sorted(expected.intersection(*(set(v) for v in estimates.values())))
        report[scene] = {"expected_epochs": len(expected), "common_epochs": len(common),
                         "methods": {}}
        for method, values in estimates.items():
            enu = np.asarray([helper.enu_error(values[s]["latitude_deg"],
                             values[s]["longitude_deg"], values[s]["height_m"], *truth)
                              for s in common])
            distance = np.linalg.norm(enu, axis=1)
            report[scene]["methods"][method] = {
                "three_d_rmse_m": float(np.sqrt(np.mean(distance**2))),
                "horizontal_rmse_m": float(np.sqrt(np.mean(np.sum(enu[:, :2]**2, axis=1)))),
                "vertical_rmse_m": float(np.sqrt(np.mean(enu[:, 2]**2))),
                "p95_m": float(np.percentile(distance, 95)),
                "gross_gt10m": float(np.mean(distance > 10)),
                "available_epochs": len(set(values) & expected),
                "missing_seconds": sorted(expected-set(values)),
            }
            for s, error in zip(common, distance):
                errors.append({"scene": scene, "method": method, "second": s,
                               "error_3d_m": float(error)})
    (OUT/"PVT_SCORING.json").write_text(json.dumps({
        "scope": "offline correction proxy; one retained seed; not closed-loop reproduction",
        "truth": truth, "scenes": report,
    }, indent=2))
    with (OUT/"PVT_EPOCH_ERRORS.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(errors[0]))
        writer.writeheader()
        writer.writerows(errors)
    lines = [
        "# GSS 现有 IQ 离线回放结果",
        "",
        "范围：一个保留噪声种子、五场景；不是独立确认实验，也不是论文闭环跟踪复现。",
        "未重新生成 IQ。扩范围 GPU 相关采用现有扫描器的副本，原接收机未修改。",
        "GSS 使用理想三角参考峰、四轮搜索、0.001 chip 容差，无新增门控。",
        "每秒窗口取 100 个相干历元估计的算术均值；不等于完整 1 秒相干积分。",
        "clean 的 27–39 秒窗口估计常数零点，40–79 秒用于评价。",
        "",
        "## 共同观测量 PVT（3D RMSE，米）",
        "",
        "| 场景 | 原接收机 | EML 0.1 chip | 仓库 DD | GSS | 共同历元 |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for scene, record in report.items():
        values = record["methods"]
        numbers = " | ".join(f'{values[m]["three_d_rmse_m"]:.4f}'
                             for m in ("original", "eml", "dd", "gss"))
        lines.append(f"| {scene} | {numbers} | {record['common_epochs']} |")
    lines += [
        "", "全部方法只修改 G09 C1C，其余观测和星历、RTKLIB 配置相同。",
        "完整可用率、缺失历元、水平/垂直 RMSE、P95 和粗差率见 PVT_SCORING.json。",
        "所有估计结果及解算输出先封存，再读取模拟器位置用于评分。",
        "",
        "## 解释范围",
        "",
        "- 保留数据已用于过往探索，单个种子不能给出跨独立运行的统计显著性。",
        "- 窄相关器和 DD 沿用仓库实现；这里比较的是统一离线修正流程。",
        "- GSS 查询值由 0.01 chip 网格线性插值得到，并非原论文按需直接相关。",
        "- 原扫描器未重放完整前端 FIR；相关器输出零点已校准，但形状误差可能仍存在。",
        "- 此处没有 DLL 闭环反馈和原论文码环滤波，不能将失败或成功外推为论文方法结论。",
        "- CORRELATOR_DIAGNOSTIC.json 的 clean 差分为同状态反事实代理，不是几何伪距真值。",
        "- 程序自检通过仅代表实现检查通过，不代表多径抑制或 PVT 获益。",
    ]
    (OUT/"REPORT.md").write_text("\n".join(lines)+"\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
