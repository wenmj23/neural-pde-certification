# 新冻结网络的两个可复现案例

这两个案例使用10月30日在规则锁定后生成的seed49网络，其中k=2是此前固定网络中没有使用过的频率。它们仍来自同一种一维Poisson方程和训练方式；目标是验证冻结结果，不重新训练。

Windows PowerShell，在仓库根目录：

```powershell
.\.venv\Scripts\python.exe -m pip install -r requirements-certifier.txt
.\.venv\Scripts\python.exe scripts\certify.py --case configs\examples\holdout_shallow_k2.json
.\.venv\Scripts\python.exe scripts\certify.py --case configs\examples\holdout_deep_k5.json
```

两个配置均使用三次残差重构、均匀网格和范围紧度目标U/L≤1.02，同时对绝对容差0.001作判断。 这次干净副本重跑：浅层k2的最大误差范围约为[0.0001836784, 0.0001845465]，判断满足；深层k5约为[0.0108158466, 0.0108224399]，判断超过。正式精确端点以输出result.json为准。网络权重文件由配置中的SHA-256固定；程序生成新的结果目录，包含有效输入、精确有理端点、网格与源码。范围紧度不是候选解相对误差。更多解释见[统一入口教程](certifier_quickstart.md)。

决定容差时可以运行：

```powershell
.\.venv\Scripts\python.exe scripts\certify.py --case configs\examples\holdout_deep_k5.json --goal decision --error-tolerance 0.01
```

连续保证针对冻结binary64权重所定义的精确数学网络，依赖模型、边界条件、导数包含和区间后端正确性，尚无独立/形式化软件审计。完整历史研究脚本的原始参考依赖与此日常入口不同。
