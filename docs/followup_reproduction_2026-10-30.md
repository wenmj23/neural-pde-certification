# 10月30日后续研究：如何查看及重跑

这次提交包含工作记录与研究脚本，并未改 `scripts/certify.py`。公共入口仍仅支持零边界、正弦源。所有范围都是依赖数学假设和区间后端的条件性保证；尚无独立/形式化审计。

## 最容易独立重跑：6个多项式源项冻结网络

以干净克隆为例，先安装现有 `requirements-certifier.txt`，不需要PyTorch重新验证。训练冻结前使用PyTorch 2.14.1+cpu，**仅重新训练**需要PyTorch；不保证其他平台上的重训权重逐位一致。

研究脚本输出到**仓库目录的同级** `run_20261030`，并拒绝覆盖已有输出。下面PowerShell命令在仓库根目录运行；请确保同级没有本实验的旧输出：

```powershell
python -m pip install -r requirements-certifier.txt
New-Item -ItemType Directory -Force ..\run_20261030\polynomial_pinn_training | Out-Null
Copy-Item configs\frozen_polynomial_pinns_2026-10-30.json ..\run_20261030\polynomial_pinn_training\frozen_polynomial_pinns.json
Copy-Item results\2026-10-30\followup\polynomial_pinn_training\fixture_sha256.txt ..\run_20261030\polynomial_pinn_training\fixture_sha256.txt
python scripts\78_certify_polynomial_pinns.py
python scripts\79_polynomial_threshold_analysis.py
```

78重新读取6个已冻结网络、构造6个已知解的独立参考与48个固定网格范围；79事后利用参考构造阈值、读取78结果，因此不算新的验证计算。新克隆中原有的18/24号参考脚本、interval_certification及依赖必须保留；78调用74的原型定义，不能删掉74或改变源码分隔符。`frequency: 1`是旧网络读取器校验的兼容占位符，不表示源项是正弦；真正方程由 `source_coefficients` 指定。已有权重SHA与训练配置列于[结果](../results/2026-10-30/followup/polynomial_pinn_training/)。

## 其他脚本的依赖

- 54–60 的部分研究依赖同级 `run_20261030/public_crosscheck/published_references.json` 和前序运行；参考可由仓库 `results/2026-10-30/new_networks/network_references.json` 复制。58输出完整组合供59/62使用；也可用66生成62所需的9个输入描述。60需要59运行时产生的源码快照。
- 61为可控多项式构造，63/64/67/68依赖61或67的原型定义；62需要58的网络输入描述或66生成的描述，76依赖62和74的定义。每个实验使用不同新输出目录。
- 69–71使用仓库原有18和24参考脚本与seed49网络；72读取71生成的参考，73读取72结果。74生成多项式源项旧网络研究结果，75依赖74的输出，76是导数数值点检。77是额外训练，78加载77的冻结输入，79分析78输出。
- 有些脚本用 `exec` 加载同仓库的研究脚本前缀，结构脆弱，属于记录当时实验的研究原型，不能当作稳定接口。需要更改时应重新验证，不要静默覆盖原实验。

[精选结果](../results/2026-10-30/followup/README.md)保留汇总、部分去除本地路径/终端信息后的逐项结果、新冻结权重和训练日志。原始历史ZIP及历史逐次快照并未公开；不同目录的统计可能使用同一网络或重复网格，不可相加视作独立样本。
