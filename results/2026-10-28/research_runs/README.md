# 2026-10-28 后续研究运行归档

此目录保留脚本 05–13 的选定运行证据，早期 01–04 结果继续保留在上一级目录。没有改写数值输出或重跑实验；原始 UTC 时间戳目录名保留。

`archive_manifest.json` 逐项记录来源生成路径、归档路径、SHA-256 和字节数。文件包括配置、原始结果、检查、单次计时、图及必要网络权重。首次固定案例加载过程的中间运行及可重复生成的较大参考曲线 CSV 未重复纳入；最终加载结果与图保留，具体选择见 omissions。

| 目录 | 内容 |
|---|---|
| `05_cross_problem_validation/` | 六次训练的权重、配置与 24 个状态 |
| `06_known_error_and_cost_checks/` | 27 个已知误差组合与早期成本 |
| `07_high_frequency_methods/` | 常数/直线、自适应、219 个配置 |
| `08_matched_tolerance_benchmark/` | 同间隙成本、容差判定 |
| `09_interval_prototype/` | 早期区间原型 |
| `10_function_generalization/` | 混合频率、边界峰等补充测试 |
| `11_frozen_main_case/` | 主案例重新加载和基准复现 |
| `12_interval_audit/` | 精确端点、连续参考和 94 项检查 |
| `13_verification_cost/` | 36 个计时请求、阶段分析、143 项检查 |

新克隆后，在项目根目录运行 `python scripts/14_restore_archived_results.py`，即可按清单恢复旧脚本依赖的 `results/generated/` 历史路径。该步骤只复制和校验，不训练或执行实验；遇到不同内容的已有文件将拒绝覆盖。大参考曲线可由脚本 11 重新生成。

早期脚本的固定生成路径可能被再次运行覆盖，归档目录不会因此变化。环境和源码哈希描述保存当时的运行，不保证跨平台耗时相同。源代码历史注释中的 working/uncommitted 及旧工作记录中的“尚未同步”指当时状态，日终归档的当前说明见根 README 与总结。
