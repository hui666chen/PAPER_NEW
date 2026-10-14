# 短时延"区间 PVT"信息量测试

协议（冻结规则）见 `protocols/SHORT_DELAY_PVT_INFORMATION_TEST_PROTOCOL.md`。

不生成 IQ、不跑 GPU 接收机，只用已有的 PRN09 短时延试点干净窗口和响应库。

## 运行

在实验电脑的仓库根目录（`PAPER_NEW`）下：

```bash
# 0) 先跑自检（不需要数据，约 2 分钟），确认环境没问题
python3 scripts/pvt_information_test/stage_a_cases.py --self-test --workers 2

# 1) Stage A：单星半合成（真实干净曲线 + 已知反射），4 个格子并行
python3 scripts/pvt_information_test/stage_a_cases.py --workers 4

# 2) Stage B：多星 PVT 蒙特卡洛（Stage A 判 STOP 时会自动拒绝运行）
python3 scripts/pvt_information_test/stage_b_pvt_mc.py

# 可选：用真实场景天空图（RTKLIB rnx2rtkp 加 -y 2 输出的 .pos.stat）
python3 scripts/pvt_information_test/stage_b_pvt_mc.py --rtklib-stat <path>.pos.stat
```

默认路径与 `refine_bias_intervals.py` 一致：
`/home/zoxoai/Documents/ChatGPT/paper/gps_l1_prn09_short_delay_pilot_20261013/calibration_v2_after_acquisition_amendment`。
路径不同就加 `--paper <paper 目录>` 或 `--data-dir <calibration 目录>`。

依赖：numpy、scipy、pandas；Stage A 会从本地导入 `receiver_chain_response_control`
（和之前的脚本一样）。

## 输出（`results/pvt_information_test/`）

| 文件 | 内容 |
|---|---|
| `STAGE_A_REPORT.md` / `stage_a_summary.json` | 区间覆盖率、宽度、单侧比例、DD/EML/中点修正误差、检测率、判定 |
| `stage_a_cases.csv` | 每个观测一行（约 1.6 万行），Stage B 的输入 |
| `STAGE_B_REPORT.md` / `stage_b_pvt_summary.json` | 各场景各方法 3D RMSE、区间法/最强对照比值及 95% CI、判定 |

把 `results/pvt_information_test/` 整个提交即可（csv 约几 MB）。

## 耗时

Stage A 每个格子约 6000 次剖面，每次约 0.06–0.15 s；4 核并行约 10–20 分钟。
Stage B 默认 4000 次 × 12 个场景，约 5–10 分钟。
