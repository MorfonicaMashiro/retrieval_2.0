# 前向 TB 空间平均算子

基线算子保留在 `collocate_day.m`：先对同一 SMOS–OIB 3×3 邻域内的
`hi/hs` 求原始点算术平均，再计算一次 TB。三种新增算子分别放在独立目录，
验证入口均使用对应的 `verify_*_operator` 名称。

正式试验统一使用 `centerStride=3`、`centerPhase=[2 2]` 的互不重叠窗口；
步长 1 结果仅作为历史记录，不再用于后续比较。

| 目录 | 算子 | 状态 |
|---|---|---|
| `baseline_parameter_mean_stride3/` | 原始参数算术均值后算 TB；中心格步长 3 | 已实现，独立输出 |
| `operator_2_sample_tb_mean/` | 逐点干舷推逐点冰厚/雪深，逐点算 TB 后算术平均并重新反演 | 已实现；旧前向脚本仅保留为诊断 |
| `operator_3_weighted_parameter_mean/` | 到 SMOS 中心的反距离加权干舷/`hi/hs` 后作单柱 TB 与反演 | 已实现 |
| `operator_4_weighted_sample_tb_mean/` | 逐点干舷推逐点冰厚/雪深，逐点算 TB 后等面积加权并重新反演 | 已实现 |

算子 3 当前使用 OIB 点到 SMOS 中心格中心的反距离权重 `w=1/d`；权重指数
可通过 `distancePower` 调整。旧的 2.5 km 等面积子网格权重仍作为
`equal_area_subgrid` 选项保留，但不再是算子 3 的默认设置。

真正逐点反演的公共数值核心放在 `common/`；每种路线的入口、结果和说明仍
分别保存在对应算子目录。所有新结果均写入各自的 `results/pointwise_inversion/`，
不覆盖 `processed/` 历史缓存。
