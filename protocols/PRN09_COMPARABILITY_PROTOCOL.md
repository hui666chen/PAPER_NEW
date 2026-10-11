# PRN9 comparable current-binary validation

冻结日期：2026-10-10（Asia/Shanghai）  
状态：运行前冻结；不得根据验证结果修改参数。

## 目的

重新验证 PRN9，使它与 PRN30 最新确认批次严格可比。旧 PRN9 结果保留为历史结果，不参与本批次参数选择，也不替代本批次。

## 固定处理链

- 接收机二进制：`<LOCAL_USER>/chh/GPU/High/gnss-sdr-cuda-scan/build-cuda-scan/src/main/gnss-sdr`。
- 校准响应库：必须使用该二进制重新生成的 PRN9 专用响应库；不得复用旧 PRN9 或 PRN30 响应库。
- GPS L1 C/A，PRN9，采样率 25 MHz，IF -4.58 MHz，401 个相关抽头。
- 共享复增益比模型；每历元增益剖面掉，反射/直达复比在窗口内共享。
- 反射延迟固定为 0.8 chip，反射相位固定为 `8.833159347728436` degree，与 PRN30 最新批次完全一致。
- 延迟估计使用与 PRN30 最新批次相同的粗网格 shared-ratio tau 平均值；不进行 PRN9 旧批次的连续局部 `±0.01 chip` 再搜索。
- 单窗口证据：`cv_gain_min > T_null` 且两个交叉验证方向均不触及候选边界。
- 持续性门控：固定 `3-of-5`，不在本批次重新选择 `K/M`。
- 阈值固定为已有独立 clean 标定批次得到的 `T_null = 0.0032589351980997073`。

## 数据和随机种子

- 校准种子：`2026101100`；仅用于当前二进制下的 PRN9 clean 响应库。
- 确认种子：`2026101110`、`2026101111`、`2026101112`、`2026101113`、`2026101114`。
- 每个确认种子各生成一次 145 s IQ；clean、0.05、0.10 三个场景从同一条 IQ/同一份跟踪状态的对应窗口提取。
- 每个场景包含 16 个 1 s 窗口：clean 为 61–76，0.05 为 100–115，0.10 为 121–136。
- PVT 阶段每个种子只运行一次源接收机，然后由同一份源 RINEX 复制出三种 PVT 观测流；只有 G09 的 C1C 字段允许不同。

## 三种 PVT 方法

1. `single`：PRN9 不加伪距修正；
2. `always_shared_ratio`：所有评分窗口使用 shared-ratio tau 平均值；
3. `gated_persistence_shared_ratio`：仅在固定 `3-of-5` 门控激活时使用该修正，否则回退到 single。

修正固定为

\[
\delta\rho_w=-\frac{c}{1.023\times10^6}
  \frac{\hat\tau_{A\to B,w}+\hat\tau_{B\to A,w}}{2}.
\]

## 评分

仅在 PVT 解算结束后使用模拟器真值评分。以 5 个确认种子为独立运行单位，报告水平、垂直、3D RMSE、P95 位置误差、10 m 粗差率、PVT 可用率以及 paired 胜/平/负。

- clean：`gated_persistence_shared_ratio` 相对 `single` 的 3D RMSE bootstrap 95% CI 下界不得低于 -0.05 m；
- 幅度比 0.05：同一非劣标准；
- 幅度比 0.10：平均 3D RMSE 改善至少 0.05 m、相对改善至少 1%，且 CI 下界大于 0。

bootstrap 只对 5 个独立确认种子的配对差值重采样，种子固定为 `2026101099`。所有中间日志、原始 GNSS-SDR 日志、响应库和哈希均保留。

## 审计锚点

- 共享规则文件：`<LOCAL_USER>/Documents/ChatGPT/paper/gps_l1_prn30_shared_ratio_v2/PERSISTENCE_RULE_FREEZE.json`，运行前 SHA-256：`aa72890c5a96ed8e424229be0b78f07771ab57bdea974cc4560ba31578b0cc3b`。
- 阈值文件：`<LOCAL_USER>/Documents/ChatGPT/paper/gps_l1_prn30_shared_ratio_v2/null_calibration/threshold_freeze.json`，运行前 SHA-256：`0126100761b160d4504d5fc63e558ae77aec96d08e072b545a5e2385db656e83`。
- 当前接收机二进制运行前 SHA-256：`1ddccdc3c7f4fde5680b38db8a139baa6d8d3692229d74e7d3d3b5e05683840b`。

真值不能参与响应校准、单窗口诊断、持续性门控或 PVT 观测生成；只允许参与最终评分。
