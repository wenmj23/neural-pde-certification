# 历史实验结果

本目录保留会话中原始已执行实验的输出，不由新运行覆盖。数值来源可追溯到 notebooks/ 中的执行输出和嵌入元数据。

- baseline_training.csv：训练检查点与参考网格误差。
- certificate_comparison.csv：固定网络、1001 残差检查点的认证对照。
- signed_results.csv：带符号方法在不同检查间距下的结果。
- final_weights.json：固定网络的训练后参数。
- console/、figures/：原始执行输出与图。

误差单位为原问题单位；ratio 为上界/密集网格参考误差。所有上界是理论公式的普通浮点评估，未作完整舍入认证。
