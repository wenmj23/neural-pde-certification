# 2026-11-02 实验结果索引

这些研究由助手在 Linux 执行，用户尚未在其 Windows 机器复现。本索引列出的数值检查受数学前提、软件及区间后端正确性约束；它们不构成形式化软件审计。详细结论和局限见[当日总结](../../reports/research_summary_2026-11-02.md)。

- [新多项式源项 PINN](new_forcing_study/)：12 个网络的训练配置、日志和冻结权重在 [training_metadata](new_forcing_study/training_metadata/)；证书、点值诊断、两种计时及存储故障对照分别位于对应子目录。权重 SHA-256 为 `28c3c7a9da832b21af635dee344cc25835e5f9fd0e74e5e83301aa0657c83795`。
- [不同频率正弦 PINN](sine_frequency_study/)：12 个网络的训练配置、日志和冻结权重在 [training_metadata](sine_frequency_study/training_metadata/)；含连续参考、训练诊断及预算失败未删的时间记录。权重 SHA-256 为 `327f29609274496cabcfb784cbd8ba1bda23a2b17b1d1b51abf82c2415030d80`。
- [正弦导数点检和人工采样反例](sine_point_audit/)：`alias/control_fixture.json` 是从一组 k=1 权重仅更改频率和标签得到的 k=7/k=77 人工控制，**不是**训练所得 PINN；`initial_precision` 保留 60/90 位精度下预期退出的不稳定数据与源代码快照。
- [多项式迁移、缩放和敏感性](extended_rounds/)、[全域导数检查](full_domain/)、[点诊断](reusable_point_audit/)、[预算判断](decision_budget/)、[同网格基线](new_pinn_baselines/)和[漏项定位](omission_localization/)：对应研究报告列出实验选择方式和局限。

从仓库根目录执行 `python scripts/145_verify_new_forcing_20261102.py --training results/2026-11-02/new_forcing_study/training_metadata --certificates results/2026-11-02/new_forcing_study/certificates`，以及 `python scripts/157_verify_sine_frequency_20261102.py --training results/2026-11-02/sine_frequency_study/training_metadata --certificates results/2026-11-02/sine_frequency_study/certificates`，可独立读回两组输入、证书键、哈希和保存的算术关系。这些检查不重跑神经网络和区间内核；更多命令及依赖见[复现说明](../../reports/reproduction_2026-11-02.md)。
