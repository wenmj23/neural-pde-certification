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

## 2026-10-28 后续查阅补充

以下为助手查阅摘要、相关方法和条件的状态；没有完成全文逐定理精读或算法完整复现。访问日期 2026-10-28。

| 原始来源 | 借鉴点与边界 |
|---|---|
| Demlow, Franz, Kopteva, Maximum norm a posteriori error estimates for convection–diffusion problems, IMA JNA 43(5), 2023, 2562–2584, https://doi.org/10.1093/imanum/drad001 | Green 函数与加权残差、自适应层现象；总误差还含额外半范数/数据振荡，条件不同，不能直接套用到本例 |
| Liu, Huang, Protopapas, Residual-based error bound for physics-informed neural networks, UAI 2023, https://proceedings.mlr.press/v216/liu23b.html | 方程结构与残差构造误差界，紧度与成本取舍；主要为多数线性 ODE、部分非线性 ODE、一阶线性 PDE，非一般椭圆边值问题的现成公式 |
| mpmath 官方 contexts 文档，https://mpmath.readthedocs.io/en/latest/contexts.html | 区间包含语义与实验性状态；当前原型锁定安装的 1.3.0 并记录源文件哈希，文档页面可能为较新开发版 |
| NIST DLMF 3.1，https://dlmf.nist.gov/3.1 | 定向舍入、区间与有理数表示的背景；不提供对本项目代码的正确性证明 |

Führer–Rojas 的残差有限维部分加剩余量框架是本轮的重要方法线索。分段直线残差加 Kh²/8 在本项目按经典插值和比较原理自行实现；尚没有首创性或完整复现某篇论文的证据。更详细的对应关系见 notes/working/high_frequency_round_2026-10-28.md。
