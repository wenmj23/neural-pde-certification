# 2026-11-01：公开冻结网络的独立目录复现

助手在 Linux、Python 3.12、python-flint 0.9.0、mpmath 1.3.0 执行；用户尚未在 Windows 复现。为减少对同级历史运行目录的依赖，新增 `scripts/106_portable_polynomial_check_20261101.py`，直接读取公开的两组冻结网络及训练配置；此脚本依赖仓库里原有的 `scripts/74_polynomial_forcing_prototype.py`、参考实现和区间后端，不需要 PyTorch 和原有结果文件。它仍为研究原型，不是公开 CLI 的多项式功能。原型使用提取 74 号脚本定义的兼容方式，如果 74 号脚本结构改变会明确报错。

每个网络在 64 个均匀区间上分别执行二次和三次证书，再独立使用 mpmath 区间参考估计连续最大误差；每个网络的参考在二次与三次证书计算完成之后计算，不参与选择。两组（种子 203/204 与 305/306）各 12 个网络，共 48 个范围；24 个参考达到指定收敛标准，48 个证书都包含对应参考区间。这是条件性数值交叉检查，仍依赖假设和区间实现正确性，不构成独立形式化审计，也不证明任意网络均有效。对浮点权重按精确实常数解释，不认证任意浮点推理。

固定 64 段时，按 `L>0 且 U/L≤1.02`，二次 18/24、三次 24/24 达标；两组各自为二次 9/12、三次 12/12。这只表示所给误差范围窄，不能解释成真实网络误差小于 2%，也不是总体速度对比。不同来源复用种子，网络不是独立同分布样本。

与此前固定 64 段的数据比对：第一组 24/24 的上下界完全一致，第二组旧记录可对照的 12/12 上界一致。另在同机隔离代码副本中重跑第一组，24 条保存记录逐字节一致；这不是跨平台或独立实现复核。两次输出均检查预期键和数量，导出时以临时文件替换并记录哈希；对已有输出目录拒绝覆盖，但不能保证所有硬件故障场景。

结果为 `results/2026-11-01/portable_reproduction/` 下的配置、明细、参考和汇总。若在仓库同步本研究脚本后复现，可在 Windows PowerShell 安装 `requirements-certifier.txt`，然后执行（输出目录应不存在）：

```powershell
.venv\Scripts\python.exe scripts\106_portable_polynomial_check_20261101.py --fixture configs\frozen_poly_new_training_2026-11-01.json --out results\generated\portable_new
.venv\Scripts\python.exe scripts\106_portable_polynomial_check_20261101.py --fixture configs\frozen_poly_holdout_training_2026-11-01.json --out results\generated\portable_holdout
```

上述 Windows 命令尚未经用户本机验证。无需联网或访问原始历史结果，但需要安装的依赖、仓库代码以及冻结网络文件。每组输出的 `summary.json` 为最终完成标志；中断时输出目录可能只包含部分结果，不能当作完整实验。需继续加强持久化方案：本脚本原子替换的是单个文件，并未提供跨文件事务或断点续跑。
