# Seed policy for the confirmatory outlier Monte Carlo

冻结时间：2026-09-25

状态：`FROZEN_BEFORE_CONFIRMATORY_MONTE_CARLO`

## 绑定的预测

预测文件：`gps_l1_prn14_blind_transfer/outlier_theory_confirmatory_prediction_27_gaussian_frozen/predictions.json`

SHA-256：`a43c47b22963bcf1ac438b08764f8a7159a56abd9ea5b51f620a46e269f67bd7`

## 每格试验次数

- 若封存的 `P_global_pred < 0.02`，该格直接运行 trial index `0..999`，总计 1000 次；
- 其余格直接运行 trial index `0..199`，总计 200 次；
- 次数只由封存预测决定，不读取前 200 次结果后再决定是否追加。

提前确定的 1000 次格点为：

- `(delay=0.70, amp=0.23, N_epoch=100)`；
- `(delay=0.90, amp=0.25, N_epoch=50)`；
- `(delay=0.85, amp=0.30, N_epoch=100)`。

## 种子派生规则

每个格点、每个 trial 使用独立确定性种子。规范字符串使用 ASCII：

```text
prn14-confirmatory-outlier-v1|
prediction_sha256=a43c47b22963bcf1ac438b08764f8a7159a56abd9ea5b51f620a46e269f67bd7|
delay={delay:.2f}|amp={amp:.2f}|epochs={N_epoch:d}|trial={trial_index:06d}
```

对完整单行字符串计算 SHA-256，取 digest 的前 16 bytes，按 4 个 little-endian `uint32` 解释并传给 `numpy.random.SeedSequence`；位生成器固定为 `numpy.random.PCG64DXSM`，再构造 `numpy.random.Generator`。

每个 trial 内随机数消费顺序固定为：

1. `normal(size=(100,401))` 生成复噪声实部；
2. `normal(size=(100,401))` 生成复噪声虚部；
3. 两者按 `(real + 1j*imag)/sqrt(2)` 合成；
4. `choice(gain_pool, size=100, replace=True)` 生成 100 个历元增益，并平方得到权重。

之后使用冻结噪声协方差的特征值/特征向量变换复白噪声。`N_epoch=20/50/100` 分别使用该 trial 生成的前 20/50/100 个历元。

## 独立性和禁止项

- 不同格点因规范字符串包含格点参数而使用不同种子；
- 同一格点不同 trial 因包含 trial index 而使用不同种子；
- 批处理大小不得影响任何 trial 的随机样本；
- 禁止试跑多个主种子后按评分结果选择；
- 禁止因中间 MC 结果修改 trial 数、格点、事件规则或预测文件；
- 原 27 格 MC 流、先前 copula 审计流和本确认实验的随机流相互独立。

冻结后不再修改本策略；任何字节变化均使对应运行清单失效。
