# 步长 3 平均次序对比

本目录对比相同 733 个非重叠 3×3 窗口上的两种前向 TB 算子：

- 参数先平均：`F(mean(hi), mean(hs))`；
- TB 先计算：`mean(F(hi_m, hs_m))`。

运行：

```matlab
cd('D:\海冰厚度联合反演\retrieval_matlab_ver\analysis\tb_spatial_averaging\comparison_stride3')
compare_stride3_averaging_results
```

`results/three_way_metric_comparison.csv` 同时列出论文报告值、历史步长 1、
步长 3 参数先平均和步长 3 TB 先计算。历史步长 1 仅作回顾，不作为严格
单因素对照。
