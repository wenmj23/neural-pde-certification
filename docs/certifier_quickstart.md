# 统一验证入口：从两个冻结案例开始

这个入口把已有区间验证内核组织成“读取案例→重构残差→计算连续误差范围→检查目标→保存证据”的流程，不重新训练网络，也不读取研究参考结果来决定网格或停止。

当前适配的问题是

\[
-u_*''=(k\pi)^2\sin(k\pi x),\quad x\in[0,1],\qquad u_*(0)=u_*(1)=0,
\]

候选解为v(x)=x(1−x)N(x)，N是冻结的全连接tanh网络。误差指标为整个区间上的最大绝对误差 \(\|v-u_*\|_\infty\)。示例覆盖一层和三层网络。

## 1. 安装与运行

推荐Python 3.12。在仓库根目录，Windows PowerShell：

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements-certifier.txt
python scripts/certify.py --case configs/examples/shallow_sine.json
python scripts/certify.py --case configs/examples/deep_sine.json
```

macOS/Linux创建虚拟环境后用 `source .venv/bin/activate`，后续命令相同。

统一入口只需python-flint==0.9.0和mpmath==1.3.0。不需要安装训练所用的PyTorch、SciPy或绘图库，也不需要恢复昨天的研究结果目录。

预编译包和支持平台说明见[python-flint官方安装页](https://python-flint.readthedocs.io/en/latest/install.html)。当前实际运行验证在Linux/Python 3.12完成，学生自己的环境仍待复现。

已归档的正式输出可直接查看：[浅层报告](../results/2026-10-29/examples/shallow/report.md)、[深层报告](../results/2026-10-29/examples/deep/report.md)。重跑时比较误差端点和判定；耗时会随机器变化。

两个案例使用已冻结的 `configs/frozen_networks_depth_2026-10-29.json`，配置中固定了该文件的SHA-256。文件改变时入口会拒绝继续，避免把不同权重的结果混在一起。

## 2. 先看懂两个输出

| 配置 | 网络/源项 | 默认方法 | 结果应怎样理解 |
|---|---|---|---|
| `shallow_sine.json` | 一层16单元，k=1，seed7 | 二次残差重构，均匀网格 | 最大误差约在[0.0001030742, 0.0001035296]，已证明低于0.001容差 |
| `deep_sine.json` | 三层、每层16单元，k=5，seed7 | 三次残差重构，均匀网格 | 最大误差约在[0.0036717448, 0.0036802636]，已证明超过0.001容差 |

表中仅为便于阅读的近似值；正式端点以result.json中的精确有理数和向外舍入小数为准。

第二个案例的否定判断也是正常完成的验证：验证器发现该近似解达不到要求。它不会通过改变网络权重来修好这个解。

先做三件事：打开配置找到频率、插值次数和容差；运行命令；打开程序打印的report.md，检查上下界与容差的位置关系。再看报告中的Z_L、Z_U和D_U怎样组合成最终范围。

## 3. 两个目标分别是什么

**`goal=enclosure`：把误差范围收紧。** 默认要求U/L≤1.02。它衡量“我们对最大误差知道得有多精确”，不表示网络解误差低于2%。两个示例默认使用此目标，并额外报告是否满足绝对容差。

**`goal=decision`：只判断是否满足容差τ。** U≤τ时判断满足；L>τ时判断超过；其余情况尚不能判定。它可能比得到很窄的范围更早完成。

```powershell
python scripts/certify.py --case configs/examples/deep_sine.json --goal decision
```

“容差判断已完成”和“范围已足够窄”可能一个成立、另一个不成立，程序分别报告。

预算用尽时，程序会保留最后一个完整有限范围；如果该范围跨越容差，就报告尚不能判定。尚未算出有限范围时，输出为空，绝不把它记为零误差。

## 4. 切换方法与预算

同一冻结网络可以切换二/三次残差重构和均匀/局部加密网格：

```powershell
python scripts/certify.py --case configs/examples/deep_sine.json --degree 2
python scripts/certify.py --case configs/examples/shallow_sine.json --degree 3 --mesh adaptive
python scripts/certify.py --case configs/examples/shallow_sine.json --gap-ratio 1.002
python scripts/certify.py --case configs/examples/shallow_sine.json --goal decision --error-tolerance 1e-4
```

本版由配置或命令明确指定方案，还没有加入自动选阶规则。三次需要更高阶导数，局部加密也有计算成本，它们并非总比二次/均匀更划算。

可以有意把预算缩小，观察“尚不能判定”：

```powershell
python scripts/certify.py --case configs/examples/shallow_sine.json --max-intervals 16 --error-tolerance 1e-4
```

该配置的范围约为[0.00007129,0.00013631]，跨越0.0001。程序保存范围并提示用尽网格预算。

常用配置项都在案例的verification对象中：

| 字段 | 含义 |
|---|---|
| degree | 2或3，残差插值次数 |
| mesh | uniform或adaptive |
| goal | enclosure或decision |
| gap_ratio | 上下界比例要求，例如字符串`"1.02"` |
| error_tolerance | 最大绝对误差容差，例如字符串`"0.001"`，也可以为空 |
| initial_intervals / max_intervals | 起始/最大区间段数 |
| max_attempts | 最大网格尝试次数 |
| time_limit_seconds | 软时间预算，例如字符串`"30"` |
| precision_bits | 区间精度，示例为103 |

阈值用字符串保存，比较时转为精确有理数，避免JSON浮点数悄悄改变边界。网格预算默认逐次翻倍，局部加密复用已有结果。时间预算在完整网格尝试之间检查，最后一轮可能超过软上限。

## 5. 保存了什么，怎样复现

默认输出为 `results/generated/certify/案例名/时间戳/`，也可用 `--out`指定一个新的目录。已有目录会被拒绝覆盖。

| 文件 | 内容 |
|---|---|
| report.md | 可读的方程、结论、过程与误差分解 |
| result.json | 精确端点、判定、预算、版本和输入/源码摘要 |
| trace.json | 每一档网格的范围、工作量和状态 |
| mesh.json | 最后一个完整范围对应的有理网格坐标 |
| case.json / fixture.json | 本次实际使用的配置和冻结权重 |
| case.original.json | 原始配置，命令覆盖之前的版本 |
| source_snapshot/ | 本次执行源码，保留历史可追溯性 |

重跑时将 `--case` 指向该结果目录中的case.json。它引用同目录的fixture.json；即使原来的研究配置路径改变，输入仍可复现。也可以使用保存的source_snapshot/certify.py作为入口，并将`--out`设为新的目录。

退出码0表示本次目标完成，可能是满足容差，也可能是已证明超过容差；退出码2表示预算内目标未完成；退出码1表示输入、依赖或执行错误。不要单靠退出码0判断解是否合格，应读取accuracy_decision。

## 6. 保证范围与研究记录

最终范围仍来自

\[
\max(0,Z_L-D_U)\le\|v-u_*\|_\infty\le Z_U+D_U.
\]

二次残差余项要求r∈C³，三次要求r∈C⁴；当前tanh网络和正弦源项满足这些条件。范围针对冻结binary64权重定义的精确数学网络，依赖数学模型、导数包围及区间后端的包含正确性，尚无独立或形式化软件审计。

本入口目前支持上述一维Poisson族和严格零边界；其他算子、区域、边界和激活会被拒绝，扩展需要另行实现并检验。

实现基于已有实验内核，没有新增“自动选择最优方案”的研究结论。方法和正负结果见[多层与高阶重构记录](../notes/working/depth_higher_order_2026-10-29.md)。

## 7. 日常复现与研究脚本的数据要求

上述两个统一入口示例只依赖仓库中的源码、配置、冻结权重和最小Python依赖。今天的完整原始研究ZIP未纳入公开提交，保留在工作区。需要历史连续参考文件的研究脚本17–32不能仅凭本次精选输出完整重放；它们的数据路径和条件见各工作记录。学生先复现两个统一入口，不需要运行这些研究脚本。

10月28日已公开的历史研究证据仍可按原说明用`python scripts/14_restore_archived_results.py`恢复，旧归档不受影响。
