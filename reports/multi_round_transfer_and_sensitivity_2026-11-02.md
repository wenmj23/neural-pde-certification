# 2026-11-02：多轮源项迁移、数值一致性与敏感性研究

研究由助手在 Linux、Python 3.12.14、python-flint 0.9.0、mpmath 1.3.0 执行；用户尚未在 Windows 重现。冻结网络、定义域、零 Dirichlet 边界和二/三次方法保持原范围。测试数据是已有 24 个网络，并非 24×3 个新训练网络。区间、参考和精确实权重的保证均以数学假设、实现和后端包含性质正确为前提；不属于形式化或独立人工审计。结果保存在 results/2026-11-02/extended_rounds/。

## 第一轮：未重训网络迁移到三个源项

将每个冻结网络分别用于 f=12x²、20x³、42x⁵，对应零边界真解分别为 x-x⁴、x-x⁵、x-x⁷；网络未针对这些目标重新训练。固定 64 段均匀网格，对每个网络和源项先算二次、三次证书，再算另一套 mpmath 区间连续参考。共 72 个网络-源项组合、144 个证书，72/72 个参考按预设停止标准收敛，144/144 证书与参考相交。原始参考中 128/144 被完整包含，另 16 个并未完整包含：保留这一事实与原始结果，不将相交直接说成参考被包住。

根据上述异常事后选出对应 8 个网络-源项组合，把参考收敛标准从相对 10^-6、绝对 10^-12 加紧至相对 10^-10、绝对 10^-15，最多 5000 次细分。8 个参考均收敛，两个次数各一份的 16/16 证书包含加紧后参考范围，且新旧参考相交。因此原始现象与参考范围偏宽相容，而不是独立证明验证器完全正确；这一步是事后针对性复核，不是预先盲测。

按 U/L≤1.02，二次 72/72、三次 72/72 达紧度，但网络没有针对新方程重训，误差可能很大。三个源项 24 个参考范围中点的中位数依次约 0.148、0.186、0.303；这些数值只是该批迁移网络与方程的参考误差，说明范围紧不代表解准确，也不能外推总体性能。同一网络在三方程上复用，样本不独立。

## 第二轮：高阶残差导数与源项导数

对上面 24 网络×3 源项，在 [0,1/32]、[7/32,1/4]、[15/32,1/2]、[31/32,1] 四个区间检查 r、r''、r'''、r''''，共 1152 组。Arb 和独立 mpmath.iv tanh 递推得到的区间全部相交；每段中点的 50/70 位 mpmath 数值微分稳定，且在两区间内。保存了两个精度下近似数的精确二进制有理数表示，便于独立读回计算包含关系；它们仍是数值近似，不是严格区间或证明。此轮四个区间不覆盖整域，原有源项的整域检查见 full_domain_and_frequency_stress_2026-11-02.md。

## 第三轮：符号和倍数缩放的变形测试

对原网络解 v 和对应 f，分别取 c=1、-1、2、-2，使 (v,f)->(cv,cf)。这些因子下最终线性层权重的 binary64 缩放被逐一检查为精确，其他层保持不变。24 个网络×2 网格（32/64）×2 重构次数×4 缩放，共 384 次。结果除以 |c| 后 384/384 与未缩放的证书上下界精确相等。这个性质可以查出部分不一致，但两边调用同一验证器，不能当成独立数值实现或对任意缩放因子的证明。

## 第四轮：容差判定边界

对迁移实验的 144 个非退化证书区间，使用精确有理数构造容差在下界以下、等于下界、内部、等于上界、上界以上的测试，再加三个零/单点边界，共 723 项。公开分类函数全部符合预期：上界≤容差判满足；下界>容差判超过；两者之间未知。容差等于上界应判满足，等于下界且上界更大时为未知。这是分类逻辑测试，不能验证证书数值是否正确。

## 第五轮：敏感性与读回

对第二轮保存的区间，事后模拟漏掉源项导数 f^(j)：由于源项是 [0,1] 上非负单项式，精确计算 f^(j) 在对应区间的范围，把它加到原残差区间，形成漏项情形的保守区间，再看是否漏掉保存的高精度中点参考。这没有执行修改后的验证器，只是对结果做保守变异。

1152 项中，288 项的相应源项导数恒为零，漏掉它们不会改变对应公式；其余 864 项中，618 项可以通过点值落在变异区间之外而发现，246 项仍被变异后的宽区间覆盖。未检出不意味着公式正确；此实验说明点值检查有盲区，不说明其余情形实际发生错误。

独立脚本 123 在进程退出且复制到精选目录后，利用可信原始输入核对 144 个证书、72 个原始参考、384 个缩放结果、16 条加紧参考比较、723 个边界结果及 1152 个导数结果的键、哈希和可复算的端点关系；脚本 126 又读回并复核 1152 个变异结果。二者不重跑或证明原始数值计算。124 号对 9 个临时副本进行受控变异，结果均符合预期；自洽的虚构宽区间重算哈希后被接受，说明这不是真实性认证。此前不完整明细的原始根因仍未确定。

## 复现

安装 requirements-certifier.txt。两组输入为 configs/frozen_poly_new_training_2026-11-01.json 和 configs/frozen_poly_holdout_training_2026-11-01.json，均使用新建输出目录。以下为 bash 示意，Windows 可调整路径与解释器：

```bash
F1=configs/frozen_poly_new_training_2026-11-01.json
F2=configs/frozen_poly_holdout_training_2026-11-01.json
python scripts/118_transferred_forcing_20261102.py --fixtures "$F1" "$F2" --out results/generated/transfer_new
python scripts/119_metamorphic_scaling_20261102.py --fixtures "$F1" "$F2" --out results/generated/scale_new
python scripts/120_tighten_transfer_reference_20261102.py --fixtures "$F1" "$F2" --input results/generated/transfer_new --out results/generated/tight_new
python scripts/121_decision_boundary_checks_20261102.py --input results/generated/transfer_new/records.json --out results/generated/decision_new
python scripts/122_transfer_derivative_crosscheck_20261102.py --fixtures "$F1" "$F2" --out results/generated/derivative_new
python scripts/125_forcing_derivative_mutation_20261102.py --input results/generated/derivative_new --out results/generated/mutation_new
python scripts/126_verify_forcing_mutation_20261102.py --input results/generated/derivative_new --mutation results/generated/mutation_new
```

123 的完整读回要求五个目录名分别为 transferred_forcing、metamorphic_scaling、tightened_reference、decision_boundaries、transferred_derivatives；已有精选结果可直接运行：

```bash
python scripts/123_verify_multiple_rounds_20261102.py --root results/2026-11-02/extended_rounds --fixtures "$F1" "$F2"
python scripts/124_multiple_round_guard_mutations_20261102.py --root results/2026-11-02/extended_rounds --fixtures "$F1" "$F2" --out results/generated/guard_new
```

本轮结果路径见当日索引。接下来优先考虑更独立的区间导数参考（例如有可靠余项的局部 Taylor 模型）、数据完整性根因排查和中期报告材料，而不是据少量正例声称普遍有效。
