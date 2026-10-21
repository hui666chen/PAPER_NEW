"""Finalize the report from an already completed weighted-correlator run."""
import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/weighted_correlator_comparison_20261021"
SCORING = json.loads((ROOT / "results/gss_replay_20261020/PVT_SCORING.json").read_text())
WEIGHTED = json.loads((OUT / "WEIGHTED_PVT.json").read_text())
CAL = json.loads((OUT / "CALIBRATION.json").read_text())
BASE_CAL = json.loads((ROOT / "results/gss_replay_20261020/PVT_CALIBRATION.json").read_text())
SCENES = ["clean", "d0p05_a0p20_p000", "d0p05_a0p20_p090",
          "d0p05_a0p20_p180", "d0p05_a0p50_p090"]

rows = []
for scene in SCENES:
    ref = SCORING["scenes"][scene]["methods"]
    original = WEIGHTED[scene]["scores"]["single"]
    methods = {
        "Original receiver": original["3d_rmse_m"],
        "Narrow EML 0.1 chip": ref["eml"]["three_d_rmse_m"],
        "Existing DD": ref["dd"]["three_d_rmse_m"],
        "Weighted E-L candidate": WEIGHTED[scene]["scores"]["gated"]["3d_rmse_m"],
    }
    for method, rmse in methods.items():
        rows.append({"scene": scene, "method": method, "3d_rmse_m": rmse,
                     "change_vs_original_pct": 100*(rmse/methods["Original receiver"]-1),
                     "common_epochs": WEIGHTED[scene]["scores"]["gated"]["common_epochs"]})

with (OUT / "COMPARISON.csv").open("w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0]))
    writer.writeheader(); writer.writerows(rows)

lines = [
    "# 现有 PRN9 相关器数据：统一 clean 校准流程后的共同观测量 PVT 比较",
    "",
    "## 实验口径",
    "",
    "- 同一批保留的 PRN9 IQ 和同一组每历元 401 抽头相关器输出；不重生成 IQ。",
    "- 所有延迟修正都用相同的 clean 窗口 27–39 标定零点；窗口 40–79 作为 PVT 评价段。每个鉴别器有自己的零点数值，避免把鉴别器的固有零偏差当成多径改正。",
    "- 每个方法只修改同一份 RINEX 的 G09 C1C，其余卫星、观测和星历不变。",
    "- PVT 使用相同 RTKLIB 配置和 40 个共同有效历元。原接收机、EML、DD 分数来自同批共同观测量回放；原接收机与本次重解算的 single 分支逐场景完全一致。",
    "- 加权 E-L 使用两组间距 0.10/0.20 chip、权重 1/−0.5；无自适应权重、无多径门控。每个 1 ms 历元取相关峰附近正斜率零点，逐秒对 100 历元取均值。",
    "- 零点：窄相关 EML %.6f chip，原 DD %.6f chip，加权 E-L %.9f chip（%.3f m）。原接收机分支不施加校正。" % (BASE_CAL["zero_chip"]["eml"], BASE_CAL["zero_chip"]["dd"], CAL["zero_chip"], CAL["zero_m"]),
    "- 反向修正量按接收机现有符号约定计算；若窗口内有任何无效历元，该秒加权 E-L 修正置零。",
    "",
    "## PVT 结果",
    "",
    "| 场景 | 原接收机 | 窄相关 EML | 原 DD | 加权 E-L 候选 | 共同历元 |",
    "|---|---:|---:|---:|---:|---:|",
]
for scene in SCENES:
    ref = SCORING["scenes"][scene]["methods"]
    vals = [WEIGHTED[scene]["scores"]["single"]["3d_rmse_m"],
            ref["eml"]["three_d_rmse_m"], ref["dd"]["three_d_rmse_m"],
            WEIGHTED[scene]["scores"]["gated"]["3d_rmse_m"]]
    lines.append(f"| {scene} | " + " | ".join(f"{x:.3f}" for x in vals) + " | 40 |")
lines += ["", "单位：米，3D RMSE。单个噪声种子，仅40个定位历元；属于探索性回放，不代表跨种子稳定结论。",
          "", "## 结果解释", "",
          "这次的加权 E-L 候选没有成为可用的改善方案：它在五个场景的 3D RMSE 都高于窄相关 EML；在 clean 和 0.5 幅度比场景也高于原接收机。幅度比 0.2、相位 0°时，它与原接收机几乎相同（差异约0.003 m）。因此不能把这个候选启用为默认多径抑制功能。",
          "",
          "现有回放中，窄相关 EML 是表现最好的对照：五个场景均优于本次加权候选，并在四个含反射场景中优于原接收机。它可以作为接收机已有多径缓解功能的工程基线；本轮没有证明新的加权算法有增益。",
          "",
          "对照文献的算法结构并不等于复现其完整 MAA-DLL。该论文主要用 BOC 信号，并且先估计反射延迟，再根据 BOC 延迟区间查表切换间距和权重；本次只比较了固定加权 E-L 的相关器级候选。",
          "",
          "本报告只回答这批现有数据上的共同观测量比较，不证明实时闭环性能、统计显著性或对其他 PRN/几何的泛化。"]
(OUT / "REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
print(OUT / "REPORT.md")
