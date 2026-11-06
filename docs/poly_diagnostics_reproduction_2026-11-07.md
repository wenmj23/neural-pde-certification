# 11 月 6 日两套多项式导数诊断：复现说明

本次精选保存两项关联实验：24 个已有冻结网络的三/四阶残差导数区间与高精度点值比较，以及常数/线性源项遗漏的高阶导数盲区。它们是有限数值诊断和人为故障对照，不是新的连续误差证书，也不是生产实现错误的报告。

## 保存了什么

- 两组 11 月 1 日冻结输入、固定配置、研究脚本 225、251–254、结果明细和原始汇总。
- `results/2026-11-05/independent_poly_split/locked_config.json` 与 `summary.json` 是 251 号脚本用来检查来源配置及研究原型哈希的两个前置文件；这**不代表** 11 月 5 日整套研究的原始明细已补齐。
- `results/2026-11-06/poly_split_diagnostic/` 含 864 条点值对照、48 条故意遗漏源项对照、锁定配置和汇总；`poly_zero_derivative_blindspot.json` 含 144 条对照，另有独立读回结果。
- 第 262 号脚本在新目录复制所需代码和配置，重跑 251–254 号并和仓库保存数据比较哈希。每个原实验脚本仍会拒绝覆盖已有的结果目录。

## 运行

依赖在 `requirements-certifier.txt` 中固定为 python-flint 0.9.0 和 mpmath 1.3.0；不需训练，也不需 PyTorch。从仓库根目录，以 Windows PowerShell 为例：

```powershell
python -m venv .venv
.venv\Scripts\python.exe -m pip install -r requirements-certifier.txt
.venv\Scripts\python.exe scripts\262_replay_poly_diagnostics_20261107.py --workdir "$env:TEMP\pde_poly_replay_20261107"
```

指定的 `--workdir` **必须不存在**；如果之前已经运行过，请换一个新目录名。脚本不覆盖历史结果，失败时保留该目录供检查。复现成功后其中 `replay_report.json` 显示三份原始文件与重跑文件的 SHA-256 是否一致；不同平台产生不同字节不自动代表数学错误，应查看实验的读回检查与运行环境。

也可在 Linux/macOS 已配置依赖后运行：

```bash
python scripts/262_replay_poly_diagnostics_20261107.py --workdir /tmp/pde_poly_replay_20261107
```

## 本次验证及边界

2026-11-07 在 Linux/Python 3.12、python-flint 0.9.0、mpmath 1.3.0 中，使用独立工作目录和单独虚拟环境重跑；依赖从本机已安装的同版本包复制到该虚拟环境，并非新联网安装；Windows 命令尚未经用户实际运行。检查结果为 864/864 点值相容、48/48 四次源项人为故障检出；常数/线性源项案例中原始残差 48/48 检出，而高阶导数 96 个检查未显示遗漏（源项高阶导数本来为零）。三个明细文件的重跑哈希与所保存原件一致。

这些冻结网络此前已研究过，重复设置不是独立新网络。高精度点值比较和独立进程读回不能证明整个区间或代码没有错误。故障是隔离测试中故意构造的，历史保存异常根因仍未查明。本说明只覆盖两项实验，不应宣称全部 11 月 5–7 日研究已完整复现。
