# neural-pde-certification

一维 Poisson 神经求解器的后验最大误差认证：学习、推导、实验与中期报告记录。

工作入口与阶段归档更新：2026-10-29。中期报告截止：2026-11-13。

先用 [统一验证入口教程](docs/certifier_quickstart.md) 复现两个冻结案例；查看 [10月29日阶段总结](reports/research_summary_2026-10-29.md)、[多层与高阶重构工作记录](notes/working/depth_higher_order_2026-10-29.md)及[当天示例与结果](results/2026-10-29/README.md)。[10月28日日终总结](reports/research_summary_2026-10-28.md)与[早期固定主案例](cases/poisson_1d_sine_k1_seed7/README.md)继续保留。

## 统一验证入口

加载冻结JSON网络，无需重新训练或恢复研究结果目录：

```bash
python -m pip install -r requirements-certifier.txt
python scripts/certify.py --case configs/examples/shallow_sine.json
python scripts/certify.py --case configs/examples/deep_sine.json
```

程序保存连续最大误差范围、容差判断、网格过程、有效输入与执行源码。支持二/三次残差重构、均匀/局部加密网格；明确区分“判断容差”和“收紧范围”两个目标。当前限定一维正弦源项Poisson、严格零边界和光滑tanh网络，保证依赖区间后端包含正确性，尚无独立或形式化审计。

## 当前进展

- 已实现二次/三次带符号残差重构、形状余项包络、连续范数区间界、均匀/局部加密网格，并支持一/二/三隐层tanh网络。
- 18个冻结网络的多路径研究轮完成7983项检查与933次达标计时请求；统一入口另完成60项流程检查。这些是有限问题上的方法、预算和重复测试。
- 两个示例分别给出最大误差约[0.0001030742,0.0001035296]和[0.0036717448,0.0036802636]，对0.001容差分别证明满足和超过。精确结果与源码见归档examples目录。
- 学生已大致理解整体思路，下一步在个人电脑复现；自动选阶、独立审阅、其他PDE及二维扩展尚未完成。

以下保留10月28日日终阶段进展，数值和性能描述属于当时版本：

- 从全域残差上界、局部导数界、空间包络，推进到保留符号的误差重构与剩余量估计。
- 使用分段直线近似残差，将局部余项改为 K h²/8，并测试正弦、混合频率、局部窄峰和多组固定 PINN 输出。
- 固定一维正弦源项、单隐层 16 神经元、seed 7 的 PINN，保存权重、训练来源、哈希与复现基准。
- 实现考虑验证计算舍入的区间版本，完成一次助手审查与可执行检查；保留精确有理数和向外十进制端点。
- 比较相同上下界紧度要求下的完整验证成本，并实现浮点预估网格与固定系数缓存的初步方案。

研究范围包括早期 `-u''=2` 基线和后续正弦源项，均在 `[0,1]` 上采用零 Dirichlet 边界。旧核心适配单隐层tanh，当前入口支持多隐层tanh；受控人工误差测试用于检验算法，不全部属于 PINN 训练实验。

## 代表结果与限制

同样 1001 个残差检查点，已知最大误差的 `A sin(20*pi*x)` 测试中，上界/真实误差约从常数残差方法的 11.18 降到直线方法的 1.16。局部加密对窄峰很有效，全域高频仍可能得到偏松的界。

固定正弦 PINN 的区间版本，1000 个区间时上界与连续真实最大误差之比约为 1.00058。该结论针对精确二进制权重定义的实函数，以本例数学条件和所用区间后端的包围正确性为前提；后端区间支持属实验性，尚无独立审查/形式化验证，也不认证任意浮点推理程序。

成本方面，先用浮点选网格再区间核验，三档要求的 wall time 中位数由约 2.22/8.42/21.74 秒降到 0.98/4.21/9.09 秒，端点保持一致。但每配置仅重复 3 次且有明显慢样本；缓存单独使用的收益不稳定。当时区间实现为秒级；10月29日版本、协议和结果见阶段总结，仍未展示比直接传统求解更便宜的优势。

所有实验由助手执行；学生理解与本地复现仍待完成。经典方法的实现和初步对照不直接作为首创理论贡献。

## 文件入口

| 路径 | 用途 |
|---|---|
| [最新阶段总结](reports/research_summary_2026-10-29.md) | 今天的方法、结论、限制与学生复现起点 |
| [最新示例结果](results/2026-10-29/README.md) | 两个完整示例、汇总图表与复现核验 |
| [上次日终总结](reports/research_summary_2026-10-28.md) | 昨天的完整结论与历史证据 |
| [研究日志](notes/daily/2026-10-28.md) | 当天工作索引；早期版本另行保留 |
| [理论笔记](notes/theory.md) | 基本界、重构余项、区间与停止规则 |
| [工作记录](notes/working/) | 高频、补充函数、区间审查、成本的详细记录 |
| [固定案例](cases/poisson_1d_sine_k1_seed7/) | 49 个权重、来源、基准和说明 |
| [文献清单](literature/reading_list.md) | 原始来源及实际查阅范围 |
| [历史结果](results/2026-10-28/) | 早期结果及本次选定运行的 JSON/CSV、图与校验清单 |
| `scripts/` | 实验、区间实现、成本测量及归档恢复 |
| `notebooks/` | 早期 4 份已执行笔记 |
| [下一步](NEXT_STEPS.md) | 当前待办，不代表已经执行 |

## 运行与复现

建议 Python 3.12、独立虚拟环境。Windows PowerShell：

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements-interval.txt
python scripts/14_restore_archived_results.py
python scripts/11_reproduce_frozen_case.py
```

macOS/Linux 激活命令为 `source .venv/bin/activate`。恢复脚本只校验并复制归档字节，恢复旧脚本需要的 `results/generated/` 历史路径，不训练、不运行实验。遇到不同内容的已有文件会拒绝覆盖。

后续可按需运行：

```bash
python scripts/07_high_frequency_methods.py
python scripts/10_function_generalization.py
python scripts/12_audit_interval_case.py
python scripts/13_verification_cost_benchmark.py
```

脚本 12、13 需要一些运行时间，计时时不要并行执行其他实验。脚本 10–13 每次创建新的时间戳目录；较早脚本 01–09 的输出路径固定，重跑前自行保留旧运行。已归档的证据不会随生成目录被覆盖。区间入口锁定本次检查的 `mpmath==1.3.0`，其他依赖与具体运行版本记录在结果中。

`requirements.txt` 是早期基础依赖，`requirements-experiments.txt` 补充 SciPy/mpmath，`requirements-interval.txt` 固定区间后端版本。新克隆无需重新训练即可加载主案例。

## 记录与同步约定

平时保存有意义的工作记录；用户在当天讨论结束时要求总结，再合并更新文档、代码和结果并同步 GitHub。每次交流不单独提交。保存推导、实际运行、计划和学生已掌握内容时分别标注，不创建后台定时任务。

本仓库公开，原始开题报告 PDF 及其中个人信息不纳入公开提交。工作脚本或历史记录中的 working/uncommitted 等措辞描述生成当时状态，日终进度以当前总结和归档索引为准。

本次阶段同步公开代码、文档与精选示例；今天的完整原始研究ZIP保留在工作区，未纳入公开提交。普通统一入口复现不依赖这些ZIP。需要历史参考输入的研究脚本与日常入口有不同的数据要求，见当天结果索引。
