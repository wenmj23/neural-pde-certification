# 2026-11-02：可复用的残差导数点值诊断

## 功能与范围

将上轮故障研究用到的逐点检查整理为 `scripts/134_point_audit_polynomial_20261102.py`。它接受公开格式的多层冻结网络 JSON，支持 `original`（各网络原始多项式源项）、`transfer`（三个固定的 12x²、20x³、42x⁵ 源项）和 `all` 三种模式。输入的源项系数必须为整数；它不是任意多项式源项的公共认证入口，也不支持任意神经网络架构或 PDE。固定点为 1/7、1/2、6/7，检查残差及二、三、四阶导数。

对每项用 Arb 计算该点区间，再用另一套 mpmath 网络表达式及数值微分分别在 60 和 90 位十进制精度下计算近似值。记录区间精确有理数端点、两个近似数的精确二进制有理数表示和一致性标志。这个诊断按 binary64 权重定义的精确实函数解释，不认证任意浮点推理执行。

## 正常运行与两种隔离负对照

助手在 Linux、Python 3.12.14、python-flint 0.9.0、mpmath 1.3.0 上执行；用户尚未在 Windows 本地复现。使用两组已有 24 个冻结网络、原源项与三个固定迁移源项：24×4×3×4=1152 个点值项。正常实现的 1152/1152 项满足精度稳定条件，且 90 位近似值落在 Arb 点区间内。

135 号脚本只在隔离测试子类中故意漏掉全部源项导数：1152 项中，792 项在对应点的源项导数非零，792/792 都按预期出现不相容；其他 360 项的对应导数恰好为零，360/360 仍相容。原有生产实现未被修改。792 次是测试项目数，不是 792 个独立网络。

额外使用 138 号脚本在一个已有网络上，把源项导数故意反号，而非漏掉；测试 48 项，其中应受影响的 30 项全部不相容，诊断函数返回非零，其他 18 项相容。这是有意构造的负对照，其 `sign_fault/fault_audit` 输出存在的 30 项 failures 是预期结果，不是正常生产运行的错误。此对照表明正常诊断路径会对该模拟错误返回非零，不能证明可以发现任意错误。

## 模式、读回与局限

在从已授权输入提取的单个冻结网络上分别试 `original`、`transfer`、`all`：正常及漏项对照分别有 12、36、48 条，三种模式的读回均通过。136 脚本在原进程结束并复制后，根据可信输入的网络标签、源项和固定点重算全部键、哈希、精确算术标志和汇总；正常与漏项对照各 1152 条均通过。136 只核对保存的数字，不重算两个数值后端。

137 号在临时副本中检查 8 种记录情况，行为均符合预期：能够发现截断、重复键、错误源项、错误参考或标志、汇总数量错误；如果人为改成极宽区间、并重新计算哈希且内部自洽，仍会被接受。这是结果一致性检查的已知局限，并非结果真实性认证。历史缺项根因仍未确认。

点值包含只能对所选的有限点与数值参考提供诊断。mpmath 数值微分不是严格区间证明；两套实现共享问题定义与权重。即使全部通过，也无法证明相邻点之间的区间导数界、完整验证器或任意网络均正确。该研究工具没有改变 `scripts/certify.py` 的公共功能，未来如要合入日常流程，应先稳定依赖、输入验证和运行成本并检查更广网络。

## 复现

安装 requirements-certifier.txt，输出目录必须不存在。示例 Windows PowerShell：

```powershell
.venv\Scripts\python.exe scripts\134_point_audit_polynomial_20261102.py --fixtures configs\frozen_poly_new_training_2026-11-01.json configs\frozen_poly_holdout_training_2026-11-01.json --out results\generated\point_normal_new --mode all
.venv\Scripts\python.exe scripts\135_point_audit_fault_injection_20261102.py --fixtures configs\frozen_poly_new_training_2026-11-01.json configs\frozen_poly_holdout_training_2026-11-01.json --out results\generated\point_fault_new --mode all
.venv\Scripts\python.exe scripts\136_verify_point_audits_20261102.py --normal results\generated\point_normal_new --fault results\generated\point_fault_new --fixtures configs\frozen_poly_new_training_2026-11-01.json configs\frozen_poly_holdout_training_2026-11-01.json
```

PowerShell 命令尚未在用户电脑验证。对应精选明细、负对照和模式测试存于 `results/2026-11-02/reusable_point_audit/`。读回路径见 [复现说明](reproduction_2026-11-02.md)。
