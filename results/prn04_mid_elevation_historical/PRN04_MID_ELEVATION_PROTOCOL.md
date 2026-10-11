# PRN4 medium-elevation independent validation protocol

状态：运行前冻结。目标 PRN4，固定时刻仰角约 49.4897°，同一 9 星几何。

## 固定内容

- GPS L1 C/A，目标 PRN4；SignalSim 时间、接收机位置、噪声、采样率和中频与 PRN30/PRN9 验证一致。
- 响应库必须由当前接收机二进制重新校准，窗口 68、108、129，401 taps。
- 共享比模型、单窗口阈值 `T_null=0.0032589351980997073` 和 3-of-5 持续性规则直接复用已冻结文件；不在 PRN4 数据上调参。
- 场景：clean、幅度比 0.05、幅度比 0.10；多径延迟 0.8 chip、相对相位 127°。
- 新 confirmation seeds：`2026101085`–`2026101089`。
- 每个场景 16 个窗口；真值只用于最终 PVT 评分。

## 三路 PVT

1. `single`：PRN4 C1C 不修正；
2. `always_shared_ratio`：所有窗口使用共享比延迟修正；
3. `gated_persistence_shared_ratio`：只在冻结的 3-of-5 规则通过时使用修正。

三路从同一份源 RINEX 复制生成，只允许 G04 C1C 不同；评分历元固定 9 颗卫星。

## 验收标准

- clean：门控相对 single 的 3D RMSE bootstrap 95% CI 下界不低于 `-0.05 m`；
- 幅度比 0.05：同一非劣标准；
- 幅度比 0.10：平均 3D RMSE 改善至少 `0.05 m`、相对改善至少 `1%`，且 CI 下界大于 0；
- 同时报告水平、垂直、3D RMSE、P95、粗差率、PVT 可用率和胜/平/负。

