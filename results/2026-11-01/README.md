# 2026-11-01 精选结果与检查

见[当日研究总结](../../reports/research_summary_2026-11-01.md)。本目录精选九组成功运行的明细和汇总，供交叉检查；含重复计时，重复不是独立网络。历史异常与修正保留在总结里，不把不完整文件当作有效结果。`integrity_audit.json` 和 `public_integrity_audit.json` 是对相同 1,356 条记录的检查输出，不是验证器的形式化证明。

- `poly_fixed_tolerance`：72 次固定容差判断。
- `matched_residual_ablation`、`shaped_remainder_ablation`：192 和 240 条同条件比较。
- `matched_acceptance_benchmark`、`holdout_benchmark`：各 288 次接受请求；`holdout_plan.json` 为训练前计划，脚本在训练后完成，并非外部预注册。
- `analytic_shape_stress`：120 条已知精确误差人工案例。
- `bump_mesh_comparison`、`sine_mesh_comparison`、`shifted_peak_comparison`：人工尖峰、振荡、偏移/多峰网格比较。
- 其他子目录仅保存所列实验的汇总，不能仅凭这些汇总复现全部实验。

独立拷贝运行 `scripts/100_analytic_shape_stress_20261101.py` 与 `scripts/103_sampling_alias_counterexample_20261101.py`，120 条范围记录精确一致，反例 6/6 通过。脚本把结果写到仓库**同级** `run_20261101`，并拒绝覆盖现有目标目录；请在隔离的临时工作区运行。依赖 `requirements-certifier.txt` 的区间库。其他 80–104 号研究脚本部分依赖同级目录的前序运行产物、引用路径或训练依赖，并不能保证新克隆直接全量重放。新冻结权重文件未随本次提交公开，因此依赖这些文件的实验不能只凭本仓库重放。94 号误报检查不作为有效工具发布，修正版本为 95 号。
