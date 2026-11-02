# 2026-11-02：全定义域导数检查、高频压力与记录保护

## 运行范围

本轮由助手在 Linux、Python 3.12.14、python-flint 0.9.0 和 mpmath 1.3.0 执行；用户未在 Windows 复现。对象为已有两组共 24 个冻结多项式源项网络；权重按精确 binary64 实常数解释。使用仓库固定版本及当前脚本，无重新训练。整段结果的可靠性仍依赖实现和数值后端，不是形式化证明或独立人工审计。

## 一、全定义域导数界

预先固定将 [0,1] 均分 32 段，对 24 个网络每段检查 r、r''、r'''、r''''，即 24×32×4=3,072 组。沿用上一轮两种实现：Arb 导数传播，以及独立编写的 mpmath.iv 指数 tanh 与 y'=z'(1-y²) 递推。3,072/3,072 组两套整段区间相交。每段中点的独立 mpmath.mp.diff 50 位和 70 位数值检查稳定，并全部落在对应两套区间中（3,072 项）。

相比上一轮的四个小区间，此次网格覆盖整个 [0,1]；但相交不是相互包含，更不是两者都正确的证明。每段只抽查中点，数值微分不是严格参考；两套实现共用网络输入与 PDE 定义。仅覆盖这些网络和问题，无法推广到所有网络或方程。记录中保存区间端点及点值检查布尔结果，未保存高精度点值本身，读回脚本不能独立重新计算它们。

## 二、人工构造的已知精确误差压力测试

设误差 e(x)=10^-4 sin(kπx)，则区间上的真实最大绝对误差精确为 10^-4，r=-e''。分别使用 k=1,4,16,32，16/32/64 段均匀网格及二次/三次重构，共 24 个结果，24/24 均包含精确误差。11/24 达到 U/L≤1.02 的范围紧度目标。它们是人工误差函数，不是 PINN，不衡量训练效果。

选择粗网格结果作为追问后，再对 k=16 和 32 用 128/256/512 段、二次/三次各重新计算，共 12 个结果，12/12 包含真实值、7/12 达到范围目标。两轮合计 36/36 包含。第二轮是根据第一轮结果追加的探索，不是最初预设的盲测。

具体例子：

| k | 重构次数 | 达到 U/L≤1.02 的最少网格数（仅在本轮测试的 16、32、64、128、256、512 中） |
|---|---:|---:|
| 16 | 二次 | 256 |
| 16 | 三次 | 128 |
| 32 | 二次 | 512 段仍未达到 |
| 32 | 三次 | 256 |

例如 k=32、512 段二次的 U/L 约 1.0719；三次约 1.0007。这说明在这些网格与构造案例上，高频导致验证范围变宽、需要更细网格；不说明三次对所有问题更快、更优，也不能从有限网格断定全局最少段数。U/L 衡量区间紧度，不是解的相对误差。

## 三、独立读回与受控故障

三组结果均在进程退出并复制到单独目录后，用 116 号脚本按预定频率、网格、次数或可信冻结输入标签核对完整键集合、数量、文件哈希、端点顺序及所保存标志，合计 3,072+24+12 条记录通过。源码哈希与各汇总记录一致。

117 号脚本在临时副本模拟 9 项：原始文件通过；截断、未更新哈希而改动、更新哈希后替换/重复网络键、颠倒端点、错误相交标志、错误汇总均被拒绝；一条数值被故意改为极宽区间并重算哈希后，因内部自洽而被接受，这是预期的已知局限。9/9 的行为符合预期，但检查器不能保证计算真实性、点值标志本身或数学正确。旧缺失的根因仍未知。

## 复现与记录

需安装 requirements-certifier.txt；输出目录必须不存在：

```bash
python scripts/113_full_domain_derivative_20261102.py --fixtures configs/frozen_poly_new_training_2026-11-01.json configs/frozen_poly_holdout_training_2026-11-01.json --out results/generated/full_domain_new
python scripts/116_verify_research_round_20261102.py full results/generated/full_domain_new --fixtures configs/frozen_poly_new_training_2026-11-01.json configs/frozen_poly_holdout_training_2026-11-01.json
python scripts/114_sine_known_error_20261102.py --out results/generated/sine_known_new
python scripts/115_sine_refinement_20261102.py --out results/generated/sine_refinement_new
python scripts/116_verify_research_round_20261102.py sine results/generated/sine_known_new
python scripts/116_verify_research_round_20261102.py refinement results/generated/sine_refinement_new
python scripts/117_result_guard_mutations_20261102.py --directory results/generated/full_domain_new --fixtures configs/frozen_poly_new_training_2026-11-01.json configs/frozen_poly_holdout_training_2026-11-01.json --out results/generated/full_guard_mutations_new
```

记录在 results/2026-11-02/full_domain、sine_known、sine_refinement 和 full_guard_mutations。
