# 短时延多径：Hatch 剩余空间 + 几何模型估计器测试

协议（冻结规则）：`protocols/TEMPORAL_SHORT_DELAY_TEST_PROTOCOL.md`

不生成 IQ。用你的 PRN09 响应库在相关器层面仿真 1 s 的 401 抽头曲线，
反射按真实几何变化（反射面高度 h，时延和相位随卫星仰角变），
仰角来自 SignalSim 场景用的 BRDC 星历。

## 运行（在 PAPER_NEW 根目录）

```bash
# 自检（合成响应、合成轨道，约 2 分钟）
python3 scripts/temporal_multipath_test/temporal_test.py --self-test --workers 2

# 正式（真实响应库 bank_39、真实星历、真实干净窗口残差），4 核约 20–40 分钟
python3 scripts/temporal_multipath_test/temporal_test.py --workers 4
```

默认路径：
- 响应库：`<paper>/gps_l1_prn09_short_delay_pilot_20261013/calibration_v2_after_acquisition_amendment/response_banks/bank_39.npz`
- 干净窗口（失配形状）：同目录下 `open_loop/<cell>/prn09_epoch_responses.csv`
- 星历：`/home/zoxoai/chh/GNSS/SignalSim-main/EphData/BRDC00IGS_R_20211700000_01D_MN.rnx`

路径不同用 `--paper`、`--bank`、`--calib-dir`、`--nav` 指定。
想先看理想情况可加 `--no-mismatch`（结果只作参考，判定以默认带失配为准）。

## 输出（`results/temporal_multipath_test/`）

- `TEMPORAL_REPORT.md`：每个场景 EML / HATCH / MODEL / MODEL_HATCH 的 RMS、h 估计误差、判定
- `temporal_summary.json`：全部数值和设置
- 判定为 `PASS_GO_TO_IQ` 时，额外生成 `signalsim_multipath_PRN*_h*_a*.json`：
  可直接放进 SignalSim 配置 `multipath` 字段的分段列表（每秒一段，时延和相位按几何变化）
