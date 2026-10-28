# neural-pde-certification

wmj graduation project — Neural PDE Solver Certification

一维 Poisson 神经求解器的后验最大误差认证：学习、实验与中期报告记录。

中期报告截止日期：2026-11-13。当前记录日期：2026-10-28。

## 当前已完成

- 手写近似函数：检查最大误差界与采样遗漏反例。
- 一个单隐藏层、硬边界 tanh PINN 的训练基线。
- 固定同一个已训练网络的 2×2 认证对照。
- 带符号残差的误差辅助问题与剩余量认证。

当前只研究 `-u''=2`、`x∈[0,1]`、`u(0)=u(1)=0`。后续正弦源项实验尚未运行。

所有实验最初由助手在本会话执行，用户仍需理解代码并在自己的环境复现。结果是普通 float64 对理论上界公式的实现，尚未完成向外舍入的机器级数值认证。经典推导与初步对照不直接作为新的理论贡献。

## 文件入口

| 目录/文件 | 用途 |
|---|---|
| `notebooks/` | 4 份包含已执行输出和图的独立实验笔记 |
| `scripts/` | 从笔记整理的可直接运行 Python 脚本 |
| `configs/baseline.json` | 基线实验配置 |
| `results/2026-10-28/` | 历史结果表、原始控制台输出、图和固定网络权重 |
| `notes/daily/2026-10-28.md` | 当天研究日志、结论和未解决问题 |
| `notes/daily/TEMPLATE.md` | 后续研究日志模板 |
| `notes/theory.md` | 核心公式、条件、证明线索 |
| `literature/reading_list.md` | 文献与阅读状态 |
| `reports/` | 报告说明；含个人信息的原始开题报告保留在下载包中 |
| `NEXT_STEPS.md` | 下一步实验和中期报告准备 |

## 运行

建议使用独立 Python 环境。Windows PowerShell 示例：

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
python scripts/02_pinn_baseline.py
python scripts/03_certificate_diagnostics.py
python scripts/04_signed_residual.py
```

macOS/Linux 环境激活命令为 `source .venv/bin/activate`。脚本自动将新图和结果写入 `results/generated/`，不会覆盖历史结果目录。

`01_manual_examples.py`、`03_certificate_diagnostics.py`、`04_signed_residual.py` 只需要 NumPy 与 Matplotlib；第 2 个脚本还需要 PyTorch。第 3、4 个脚本嵌入原始网络权重，因此不依赖重新训练结果。

如果使用 Jupyter，可额外安装 `jupyterlab` 并打开 `notebooks/`。运行笔记时，图和输出显示在笔记中。

## 每次有意义进展后的记录规则

1. 保存对应代码、配置、结果和图。
2. 写入当天日志：做了什么、为什么做、结论、局限、下一步。
3. 保留过去结果，新的实验使用新的目录或编号。
4. 分开记录：计划、已经运行、已经理解、尚未验证。
5. 将相关变更作为一次说明清楚的 GitHub 提交保存。

这是研究过程中的记录流程，不是后台定时任务。

## 当前最重要的结果

固定同一个网络、同样 1001 个残差检查点，基线误差界相对密集网格参考误差约为 52.49 倍；局部导数界加空间护栏约为 5.87 倍；带符号误差估计加剩余认证约为 1.034 倍。完整数值见结果表。

参考误差来自已知精确解在密集网格上的比较，不应直接称为精确连续最大误差。认证公式本身不使用精确解。

## GitHub

仓库：https://github.com/wenmj23/neural-pde-certification 。本仓库公开；原始开题报告 PDF 保留在个人下载包，不随此次提交上传。
