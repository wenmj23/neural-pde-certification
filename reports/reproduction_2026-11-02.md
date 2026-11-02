# 2026-11-02 结果读回与复现

当前文档列出新克隆仓库后从仓库根目录可运行的读回命令。实验由助手在 Linux/Python 3.12.14 执行；这些命令尚未在用户 Windows 电脑验证。读回检查数量、键、哈希和记录中可复算的关系，**不重新证明或审计数值后端**。可直接用 Python 标准库执行所列读回脚本；要重新计算证书需 `pip install -r requirements-certifier.txt`，重新训练还需 `requirements-experiments.txt`。

## 两组新训练网络与证书

以下路径均含原始冻结网络、计划配置、日志及摘要；修改网络、配置或文件格式会造成 SHA-256 不匹配：

```bash
python scripts/145_verify_new_forcing_20261102.py --training results/2026-11-02/new_forcing_study/training_metadata --certificates results/2026-11-02/new_forcing_study/certificates
python scripts/157_verify_sine_frequency_20261102.py --training results/2026-11-02/sine_frequency_study/training_metadata --certificates results/2026-11-02/sine_frequency_study/certificates
python scripts/147_verify_new_forcing_timing_20261102.py --training results/2026-11-02/new_forcing_study/training_metadata --certificates results/2026-11-02/new_forcing_study/certificates --timing results/2026-11-02/new_forcing_study/fixed_mesh_timing
python scripts/149_verify_time_to_goal_20261102.py --training results/2026-11-02/new_forcing_study/training_metadata --certificates results/2026-11-02/new_forcing_study/certificates --timing results/2026-11-02/new_forcing_study/time_to_goal
python scripts/160_verify_sine_diagnostics_20261102.py --training results/2026-11-02/sine_frequency_study/training_metadata --certificates results/2026-11-02/sine_frequency_study/certificates --grid results/2026-11-02/sine_frequency_study/grid_diagnostics --replay results/2026-11-02/sine_frequency_study/adam_replay
python scripts/162_verify_sine_time_to_goal_20261102.py --training results/2026-11-02/sine_frequency_study/training_metadata --certificates results/2026-11-02/sine_frequency_study/certificates --timing results/2026-11-02/sine_frequency_study/time_to_goal
```

重新计算例如 `python scripts/144_new_forcing_certify_20261102.py --training results/2026-11-02/new_forcing_study/training_metadata --out results/generated/poly_recompute`；生成目录必须事先不存在。正弦对应 `scripts/156_sine_frequency_certify_20261102.py` 和其 training_metadata 目录。不同运行环境可能影响耗时，不把本地重训视为逐位相同的权重。

## 点值、反例和方法对照

```bash
python scripts/136_verify_point_audits_20261102.py --normal results/2026-11-02/new_forcing_study/point_audit --fault results/2026-11-02/new_forcing_study/point_fault --fixtures results/2026-11-02/new_forcing_study/training_metadata/frozen_networks.json
python scripts/152_verify_matched_baselines_20261102.py --training results/2026-11-02/new_forcing_study/training_metadata --certificates results/2026-11-02/new_forcing_study/certificates --directory results/2026-11-02/new_pinn_baselines/accuracy
python scripts/154_verify_baseline_timing_20261102.py --training results/2026-11-02/new_forcing_study/training_metadata --baseline results/2026-11-02/new_pinn_baselines/accuracy --timing results/2026-11-02/new_pinn_baselines/timing
python scripts/166_verify_sine_point_audits_20261102.py --training results/2026-11-02/sine_frequency_study/training_metadata --normal results/2026-11-02/sine_point_audit/normal --fault results/2026-11-02/sine_point_audit/omission --alias results/2026-11-02/sine_point_audit/alias
python scripts/167_precision_incident_readback_20261102.py --initial results/2026-11-02/sine_point_audit/initial_precision/normal --recheck results/2026-11-02/sine_point_audit/alias/normal --out results/generated/precision_readback_new
```

`alias` 中的 k=7/k=77 输入是重用 k=1 权重、仅改变频率和标签的**人工控制**，不是针对这两个频率训练的网络。`initial_precision` 保存了初始 60/90 位稳定性失败，不能混入正常通过结果。

## 前期冻结网络的迁移与预算

```bash
python scripts/123_verify_multiple_rounds_20261102.py --root results/2026-11-02/extended_rounds --fixtures configs/frozen_poly_new_training_2026-11-01.json configs/frozen_poly_holdout_training_2026-11-01.json
python scripts/126_verify_forcing_mutation_20261102.py --input results/2026-11-02/extended_rounds/transferred_derivatives --mutation results/2026-11-02/extended_rounds/forcing_mutation
python scripts/140_verify_transfer_budget_20261102.py --fixtures configs/frozen_poly_new_training_2026-11-01.json configs/frozen_poly_holdout_training_2026-11-01.json --references results/2026-11-02/extended_rounds/transferred_forcing --directory results/2026-11-02/decision_budget/fixed
python scripts/142_verify_near_threshold_20261102.py --input results/2026-11-02/decision_budget/fixed --directory results/2026-11-02/decision_budget/near_threshold
```

所有报告中依赖本地历史临时目录的更早脚本可能仍不能在新克隆完整重跑；以上仅覆盖列出的本次精选数据。近门槛阈值根据已有参考事后构造，不能作为盲测。
