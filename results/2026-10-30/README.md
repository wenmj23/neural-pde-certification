# 2026-10-30：精选结果和可追溯说明

这些子目录保存本日各轮实验的配置、性能摘要、验证统计，以及两次不同次数的决策结果和参考范围。每个子目录的 `config.json` 包含运行时版本、预算和输入/源码摘要；`summary.json` 是每个案例各策略的三次计时摘要，`additional_audit.json`/`combined_audit.json` 是事后检查统计。

- [网格与增长研究](mesh_growth/summary.json)、[嵌套增长](nested_growth/summary.json)、[检查频率](check_stride/summary.json)。
- [预测策略开发](forecast_development/summary.json)和[补充评价](forecast_evaluation/summary.json)。
- [新网络连续参考](new_networks/network_references.json)、[新案例对照](new_cases_seven_methods/summary.json)、[二次三次与网格比较](degree_mesh_comparison/summary.json)。
- [普通容差对照](decision_degrees/records.json)、[误差附近约1%的门槛测试](near_threshold_one_percent/effective_settings_audit.json)、[约0.01%和256段预算测试](near_threshold_tight_budget/audit.json)。近门槛是参考引导的事后压力测试，并非盲测。

复制的 `config.json` 和 `summary.json` 来自对应正式运行，未改写原数据；为避免无效路径和冗长终端输出，决策 `records.json` 去掉了指向未公开生成目录的 `output_dir` 和 `stdout`。完整逐次试跑、网格、源码快照与历史原始归档仍保留在本地，未完整公开；这些精选数据不能替代完整逐次数据，也不足以直接运行所有历史研究脚本。`near_threshold_one_percent` 保留当时实际执行源码，因原记录 `cap=null` 的标注错误而另存有效设置审查，实际16次均为4096段上限。

日常入口可以直接运行两个[新seed49网络案例](../../docs/holdout_cases_2026-10-30.md)，不依赖历史参考目录或PyTorch。研究脚本40–53部分依赖未公开的 `results/generated/` 旧路径和原始参考，不能承诺新克隆可完整重跑研究；47号生成脚本还会拒绝覆盖已有冻结文件。范围是连续区间条件性保证，依赖光滑性和区间后端正确性，尚无独立或形式化审计。

[本次同步前的干净副本核验](sync_validation.json)：两个旧案例精确端点和源码摘要一致，两个新案例覆盖各自独立参考；新增容差拒绝、预算不足路径正常；44项选阶和100项网格流程检查通过。此前发布的208个文件按Git对象摘要逐一核对。

## 新窗口的后续成果

本日后续的组合检查、精确可控测试、边界与多项式源项研究、6个新冻结网络与训练日志，见[后续精选结果](followup/README.md)及[研究总结](../../reports/research_summary_2026-10-30_followup.md)。没有把研究原型并入公共入口，也没有上传完整历史ZIP。
