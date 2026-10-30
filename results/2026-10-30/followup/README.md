# 10月30日后续研究：精选数据

本目录是新窗口后续工作的精选数据，补充上一轮 [10月30日索引](../README.md)。

- 公共入口组合与预算：`public_crosscheck/`、`frozen_portfolio/`、`precision_crosscheck/`、`portable_replay/`、`budget_curves/`、`budget_followup/`。
- 人工可控测试与可靠性自查：`polynomial_audit/`、`derivative_crosscheck/`、`remainder_shape_check/`、`guard_checks/`、`exact_shape_identities/`。
- 非零或不确定边界研究：`boundary_lift/`、`uncertain_boundary/`、`network_boundary_lift/`、`boundary_network_comparison/`、`network_uncertain_boundary/`、`boundary_interval_envelope/`、`boundary_envelope_timing/`。
- 多项式源项和新训练网络：`polynomial_forcing/`、`poly_reference_followup/`、`polynomial_derivative_checks/`、`polynomial_pinn_training/`、`polynomial_pinn_certificates/`、`polynomial_threshold_analysis/`。[冻结网络](../../../configs/frozen_polynomial_pinns_2026-10-30.json)，格式和SHA供重放，实际源项由系数指定。

`summary.json`为正式运行汇总，部分 `records.json` 发布选定逐项数据；为移除本地绝对路径和终端输出，精简 `stdout`、`stderr`及运行路径字段，其余可读数据保留。不是完整原始归档或所有中间试跑/源码快照。近门槛根据参考事后构造；参考仅作检查；模型、点检、重复计时及人工案例不能冒充独立网络。起初2个多项式案例未完整包含较宽参考，紧参考追查见 `poly_reference_followup/`。首次导数脚本TypeError在修复后重跑，记录见[后续总结](../../../reports/research_summary_2026-10-30_followup.md)。

按[复现说明](../../../docs/followup_reproduction_2026-10-30.md)操作，部分研究依赖前序输出和专用位置。之前自动审批拒绝公开完整历史ZIP；本次未上传该ZIP、旧失败/执行快照或个人开题报告。

[同步前核验](sync_validation.json)：固定权重散列一致，21份当时记录了脚本散列的摘要与发布脚本一致；第57号实验当时未记录脚本散列。干净代码副本在同一环境重跑48个多项式范围逐项一致（不是全新Git克隆或跨平台复现）。
