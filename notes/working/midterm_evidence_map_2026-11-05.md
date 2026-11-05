# 中期报告证据索引与使用范围

> **2026-11-05 精选同步范围：**本文件与今天的研究摘要、报告草稿及接续说明已同步；只另公开 224 号数据核对脚本和其一份精简结果。11 月 5 日其余原始实验数据、配置、代码快照和 214–223 号脚本未随此次上传；下文涉及这些文件的路径是本地追溯线索，不能声称仅凭当前 GitHub 新克隆可完整复现。文中“尚未同步”等表述记录当时状态，以本说明为准。

## 可选的核心结果

| 报告主张 | 对应记录 | 使用范围 |
| --- | --- | --- |
| 一维 Poisson 残差可通过 Green 核得到 1/8 全局误差上界 | `notes/working/midterm_math_review_2026-11-05.md`；`scripts/224_midterm_evidence_check_20261105.py` | 数学推导依赖零边界和可靠全域残差界；若干有理点是辅助检查，不是证明本身 |
| 代表浅层网络可给出满足、超过、未知三种结果 | `results/2026-11-05/public_cli_replay/`；`notes/daily/2026-11-05-reproduction.md` | 今日 Linux 重跑，使用一个已公开冻结网络；用户曾在 Windows 跑较早案例 |
| 同网格方法比较的 64 段上界/参考中点中位数约 149.56、20.18、1.019、1.0003 | `results/2026-11-02/new_pinn_baselines/accuracy/records.json`；`results/2026-11-05/midterm_evidence_check.json` | 使用 12 个既有多项式网络和具体基线；中点为近似参考，基线 0/12 未证明满足不是证明超标 |
| 形状余项相比常数余项同网格 72 对更窄 | `results/2026-11-05/shape_budget/`；`reports/research_rounds_2026-11-05.md` | 同一实现中取两种上界较小值；不是独立盲测，不宣称所有情形严格更好 |
| 人工 k=32 高频案例需要更多 PDE 段数收紧范围 | `results/2026-11-05/sine_mesh_followup.json`、`sine32_budget_followup.json` | 解析构造非 PINN；1024 段以上事后选择，说明现象而非普遍性能结论 |
| 新的报告表格与数学代数检查 | `results/2026-11-05/midterm_evidence_check.json` | 保存数值和精确代数的复核，不是重新跑数值内核、也不是形式化审计 |

11 月 3 日及 11 月 5 日部分报告和明细目前只在本地工作区；GitHub 截至精选同步只包含 11 月 3 日摘要，并不完整包含这些数据。用户要求不要将历史大提交全部上传。

## 2026-11-05 核对的原始论文元数据

1. Eiras et al., **Efficient Error Certification for Physics-Informed Neural Networks**, ICML 2024, PMLR 235, pp. 12318–12347, <https://proceedings.mlr.press/v235/eiras24a.html>。摘要明确研究连续域的 PINN 残差界；不能把残差界直接等同于本项目的解误差范数界。
2. Hillebrecht and Unger, **Certified machine learning: Rigorous a posteriori error bounds for PDE defined PINNs**, arXiv:2210.03426, 2022；arXiv 页面另列 IEEE TNNLS 2023 期刊记录，<https://arxiv.org/abs/2210.03426>。此处只据摘要概述后验上界及应用方程，不据摘要断言其具体定理与本项目一致。
3. Berrone, Canuto and Pintore, **Solving PDEs by Variational Physics-Informed Neural Networks: an a posteriori error analysis**, 2022, <https://arxiv.org/abs/2205.00786>。摘要说明其对象是椭圆边值 VPINN 的能量范数估计，与本项目的 L∞ 目标不同。
4. Yu, Ho and Wang, **A Conformal Prediction Framework for Uncertainty Quantification in Physics-Informed Neural Networks**, arXiv:2509.13717v1, 2025, <https://arxiv.org/abs/2509.13717v1>。它研究统计预测覆盖；本项目简化模拟不等于论文复现，边际覆盖不能直接视作全域最坏误差证书。
5. Raissi, Perdikaris and Karniadakis, **Physics Informed Deep Learning (Part I): Data-driven Solutions of Nonlinear Partial Differential Equations**, arXiv:1711.10561, 2017, <https://arxiv.org/abs/1711.10561>。此处引用的是明确标注的预印本，不冒充核验过其他版本的完整书目信息。

以上核对的是原始文献页面标题、作者、摘要与出版元数据；不代表已逐篇逐条验证全部证明。正式提交前按学院模板调整参考文献格式，并结合所引用的具体段落核对。
