# 文献记录

状态说明：以下是助手在本会话的检索/核对状态，不表示学生已精读。

| 文献 | 原始链接 | 当前核对与适用范围 |
|---|---|---|
| Berrone, Canuto, Pintore, Solving PDEs by VPINNs: an a posteriori error analysis, 2022 | https://arxiv.org/abs/2205.00786 | 核对摘要及相关正文；椭圆问题、能量误差、残差/损失/数据项；不是直接的最大误差对照 |
| Eiras et al., Efficient Error Certification for PINNs, ICML 2024 | https://proceedings.mlr.press/v235/eiras24a.html | 核对原论文；网络导数与全域残差界，不应自动等同于解误差界 |
| Hillebrecht, Unger, Rigorous a Posteriori Error Bounds for PDE-Defined PINNs | https://arxiv.org/abs/2210.03426 | 核对摘要和相关内容；区别于二位作者早期 ODE 工作 |
| Cao et al., Residual-based error correction for neural operator accelerated infinite-dimensional Bayesian inverse problems, JCP 2023 | https://arxiv.org/abs/2210.03008 | 核对正文相关段落；借鉴误差修正思想，不是当前实验的完整复现 |
| Ernst, Rekatsinas, Urban, A posteriori certification of PDE approximations with particular application to neural networks, 2025/2026 | https://arxiv.org/abs/2502.20336 | 核对 v6（2026-08-28）相关正文；自然范数、Riesz 辅助问题、上下界 |
| Fanaskov et al., Astral: training PINNs with error majorants, 2024/2026 | https://arxiv.org/abs/2406.02645 | 核对 v2；ICLR 2026 AI&PDE workshop；能量范数实验倍数不能直接与当前 L∞ 比较 |
| Führer, Rojas, A posteriori analysis of neural network approximations, 2025 preprint | https://arxiv.org/abs/2507.06017 | 核对相关正文；离散投影残差加剩余估计，关注求积假设和常数 |
| Anwer et al., Pointwise Error Estimates for Numerical PINNs, 2026 preprint | https://arxiv.org/abs/2607.03431 | 仅核对检索返回的原始摘要；全文获取失败，定理待进一步核查 |

读每篇论文时记录：PDE 类别、求解器、误差范数、认证公式、假设、常数计算、数值求积/舍入处理、计算成本、局限。当前还没有完整复现上述论文。
