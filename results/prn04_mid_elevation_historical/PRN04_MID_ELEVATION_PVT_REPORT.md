# PRN4 中仰角独立验证报告

## 1. 目的

在不修改当前接收机、响应模型、存在性阈值或持续性规则的前提下，选择约 (40^\circ\!\text{--}50^\circ) 仰角卫星，验证冻结的共享比模型与 3-of-5 门控是否能安全服务于 PVT。

目标卫星为 GPS L1 C/A PRN4，固定几何下仰角为 **49.4897°**。采用 5 个独立噪声种子 `2026101085`–`2026101089`，每个场景 16 个评分历元。

## 2. 冻结设置

- 响应库：针对当前 GNSS-SDR 二进制重新校准，窗口 68、108、129，401 taps。
- GPU 前向形状误差：(3.3276\times10^{-5})、(3.4065\times10^{-5})、(4.0490\times10^{-5})，均通过 1% 检查。
- 反射场景：clean、幅度比 0.05、幅度比 0.10；延迟 0.8 chip，相对相位 127°。
- 单窗口阈值：复用已冻结的 `T_null=0.0032589351980997073`。
- 持续性规则：复用已冻结的 **3-of-5**，未在 PRN4 数据上调参。
- 三路 PVT：`single`、`always_shared_ratio`、`gated_persistence_shared_ratio`。
- 三路均由同一份源 RINEX 生成；只允许目标卫星 G04 的 C1C 修正不同。

## 3. 门控结果

5 个种子 × 3 个场景 × 16 个窗口，共 240 个窗口中：

- `single_window_exists` 通过数：0/240；
- 3-of-5 `gate_active`：0/240；
- clean、0.05、0.10 三个场景的门控激活率均为 0%。

原因是 PRN4 的共享比交叉验证解稳定落在分离参数边界，`cv_no_boundary=false`，因此冻结规则安全拒绝启用双径修正。

## 4. PVT 汇总

下表是 5 个种子的等权平均；RMSE、P95 和粗差率均按每个种子的 16 个评分历元计算后再平均。

| 场景 | 方法 | 水平 RMSE (m) | 垂直 RMSE (m) | 3D RMSE (m) | P95 3D (m) | 粗差率 | PVT 可用率 |
|---|---|---:|---:|---:|---:|---:|---:|
| clean | single | 2.235 | 3.055 | 3.797 | 6.175 | 0.000 | 100% |
| clean | always shared | 3.354 | 4.670 | 5.755 | 9.342 | 1.25% | 100% |
| clean | gated | 2.235 | 3.055 | 3.797 | 6.175 | 0.000 | 100% |
| 0.05 | single | 2.660 | 2.928 | 3.980 | 6.557 | 1.25% | 100% |
| 0.05 | always shared | 3.709 | 4.760 | 6.044 | 9.742 | 8.75% | 100% |
| 0.05 | gated | 2.660 | 2.928 | 3.980 | 6.557 | 1.25% | 100% |
| 0.10 | single | 2.284 | 3.410 | 4.108 | 6.425 | 0.000 | 100% |
| 0.10 | always shared | 4.351 | 6.367 | 7.716 | 11.526 | 18.93% | 96.25% |
| 0.10 | gated | 2.284 | 3.410 | 4.108 | 6.425 | 0.000 | 100% |

相对 `single` 的配对 3D RMSE 差值（正值代表改善）：

| 场景 | gated vs single | 95% bootstrap CI | 胜/平/负 |
|---|---:|---:|---:|
| clean | 0.000 m | [0.000, 0.000] m | 0/5/0 |
| 0.05 | 0.000 m | [0.000, 0.000] m | 0/5/0 |
| 0.10 | 0.000 m | [0.000, 0.000] m | 0/5/0 |

作为反事实对照，`always_shared` 的 3D RMSE 相对 `single` 分别恶化 1.958 m、2.064 m、3.608 m，5 个种子均为负向；这说明在这颗 PRN4 上强行启用双径修正是不安全的，而冻结门控成功阻止了该错误修正。

## 5. 结论范围

PRN4 的中仰角验证支持以下结论：在当前响应库和固定 3-of-5 规则下，方法在 49.49° 几何处能够稳定回退到单径，clean、0.05 和 0.10 场景的 PVT 均未劣于单径基线。

本批次**没有证明中仰角下的 PVT 精度提升**，因为门控没有激活；它证明的是跨仰角的安全性/风险控制，而不是所有仰角都能获得正向修正收益。后续论文中应将 PRN4 定位为“中仰角安全回退验证”，不要写成“中仰角改善验证”。

## 6. 审计文件

- 响应校准：[CALIBRATION_AUDIT_25M.json](CALIBRATION_AUDIT_25M.json)
- 中仰角协议：[PRN04_MID_ELEVATION_PROTOCOL.md](PRN04_MID_ELEVATION_PROTOCOL.md)
- 门控确认：[confirmation_observation/CONFIRMATION_SUMMARY.json](confirmation_observation/CONFIRMATION_SUMMARY.json)
- PVT 确认：[persistence_pvt_confirmation/PVT_CONFIRMATION_SUMMARY.json](persistence_pvt_confirmation/PVT_CONFIRMATION_SUMMARY.json)
- PVT 运行脚本：[run_prn04_persistence_pvt.py](run_prn04_persistence_pvt.py)
