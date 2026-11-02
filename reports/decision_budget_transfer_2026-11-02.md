# 2026-11-02：迁移源项的预算与容差判断

## 设计和范围

由助手在 Linux、Python 3.12.14、python-flint 0.9.0 和 mpmath 1.3.0 环境运行；用户尚未在 Windows 复现。对先前 24 个冻结网络分别使用三类多项式源项（12x²、20x³、42x⁵），共 72 个网络-源项组合；网络没有为新方程重新训练。分别在固定的 8、16、32、64 段均匀网格运行二次和三次验证器，共 576 个证书。二/三次及不同网格复用同一问题，不能当作独立网络。运行后用之前 118 号实验计算的、逐案匹配的 mpmath 区间参考核对，参考不参与证书计算。

容差在已有参考结果可见后选择，本实验不是盲测或正式预注册。参考、证书均依赖各自数值后端和假设；逐项相容不是严格软件审计。不是公开 CLI 的多项式功能。

## 固定绝对容差

对 0.1、0.2、0.3 三个容差，每个证书各判一次，合计 1728 次判断。每种次数和网格组合的 72 个结果如下：

| 绝对容差 | 证明满足 | 证明超过 | 未能判断 |
|---:|---:|---:|---:|
| 0.1 | 16 | 56 | 0 |
| 0.2 | 40 | 32 | 0 |
| 0.3 | 48 | 24 | 0 |

这些数量在本实验的每一种次数、网格下相同，1728 个判断均与参考判断兼容。阈值离这些案例的误差较远，因而 8 段就能作出判断；这不能说明低预算在其他或近门槛问题上足够。容差为绝对误差阈值，不是相对误差比例。

## 近门槛压力测试

用先前参考区间中点 m 事后构造每个案例的容差 m(1±δ)，δ 分别为 0.1%、0.01%、0.001%；对已有 576 个证书重新做 3456 次判断，没有重新运行证书。该阈值利用了参考信息，仅用于诊断最难判的门槛附近情况，不能当成实际未知真解时可自动获得的测试，也不是盲测。

以 δ=0.01% 为例，每格有 72 个组合×两个方向=144 次判断：

| 网格段数 | 二次未能判断 | 三次未能判断 |
|---:|---:|---:|
| 8 | 134/144 | 128/144 |
| 16 | 128/144 | 96/144 |
| 32 | 80/144 | 40/144 |
| 64 | 8/144 | 0/144 |

参考区间对全部 3456 个近门槛阈值均可判断；证书有判断时，与参考没有冲突。δ=0.001% 更接近误差时，64 段仍有二次 102/144、三次 56/144 未能判断。两种次数均返回“未知”是信息不足，并非证明解不满足容差。上述现象仅说明本批数据与参数的结果，不能推断三次在所有问题上更好或更快；本实验没有进行计时比较。

## 记录复核和复现

进程退出后复制独立目录，再用 140 号脚本从可信冻结输入和 118 号参考重算 576 个组合键、哈希、1728 个判断及汇总，576 个范围均与参考相交；其中 64 段的 144 个端点与 118 号此前记录逐项完全一致。142 号脚本按前述证书精确重算 3456 个事后阈值、判断和分组统计，参考未判定和不相容均为零。此流程检查保存的数据及算术关系，不重新证明区间计算正确，也不验证结果的真实性。历史缺失明细的根因依然未知。

安装 requirements-certifier.txt，使用新输出目录：

```bash
python scripts/139_transfer_decision_budget_20261102.py --fixtures configs/frozen_poly_new_training_2026-11-01.json configs/frozen_poly_holdout_training_2026-11-01.json --references results/2026-11-02/extended_rounds/transferred_forcing --out results/generated/transfer_budget_new
python scripts/140_verify_transfer_budget_20261102.py --fixtures configs/frozen_poly_new_training_2026-11-01.json configs/frozen_poly_holdout_training_2026-11-01.json --references results/2026-11-02/extended_rounds/transferred_forcing --directory results/generated/transfer_budget_new
python scripts/141_near_threshold_stress_20261102.py --input results/generated/transfer_budget_new --out results/generated/near_threshold_new
python scripts/142_verify_near_threshold_20261102.py --input results/generated/transfer_budget_new --directory results/generated/near_threshold_new
```

结果保存在 results/2026-11-02/decision_budget/。
