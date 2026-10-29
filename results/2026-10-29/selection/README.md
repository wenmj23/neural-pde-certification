# 10月29日晚：选阶、网格与学生复现精选证据

- [自动选阶汇总](degree_digest.json)：33个对象、1116次达标请求；包含开发过程重复，非独立样本数。未显示稳定优于固定策略。
- [锁定规则后6个新网络的汇总](degree_fresh_summary.json)：各策略3次重复及成本对照。
- [网格研究汇总](mesh_digest.json)：40个对象、1080次达标请求；5523个有限中间范围核对、60种最终网格与原内核端点一致。
- [代表网格对照](mesh_representative_comparisons.json)：同一阶段、同一目标的中位数与最小/最大耗时。
- [流程检查](procedure_checks.json)：选阶44项、网格100项及当时源码摘要；不等于形式化审计。
- [学生本地复现](student_reproduction.json)：依据学生贴出的终端结果记录，区别于助手Linux计时。

精选文件从已完成结果逐项复制/筛选，未重新测量或合并不同批次倍数。完整历史试跑、失败记录和源码快照继续本地保存，未公开完整归档。完整来源与适用条件见[选阶笔记](../../../notes/working/degree_selection_2026-10-29.md)和[网格笔记](../../../notes/working/mesh_selection_2026-10-29.md)。

新克隆可运行 `python scripts/34_degree_selection_checks.py` 和 `python scripts/37_mesh_selection_checks.py`（先安装 requirements-native.txt）。完整批量研究脚本33、35、36、38、39还依赖未公开的历史参考目录，不提供缺失数据即可完整复现的承诺；日常 `scripts/certify.py` 不依赖它们。优先队列和自动选择尚未加入日常入口。

本次上传前另在不含历史生成数据的干净发布副本完成[核验](sync_validation.json)：44项选阶、100项网格检查通过，浅层/深层完整误差分解的精确有理端点及公共源码摘要与原示例一致。
