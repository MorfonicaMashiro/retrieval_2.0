# 步长 3 的原始参数均值基线

本目录只改变 SMOS 中心格抽样步长：保留行列相位均为 2 的中心格
`(2,2), (2,5), ...`，因此相邻中心使用互不重叠的 3×3 EASE 窗口。

前向与反演口径保持当前基线：

- 窗口内先对 OIB `hi/hs` 作原始点算术平均，再计算一次 TB；
- Burke 递推；
- `surfaceTempC = -30 °C`；
- 不启用 Vb 进入经验介电公式前的 70 截断；
- 不增加 SMOS 质量阈值，也不按 TB 残差删点。

全量运行：

```matlab
cd('D:\海冰厚度联合反演\retrieval_matlab_ver\analysis\tb_spatial_averaging\baseline_parameter_mean_stride3')
run_stride3_parameter_mean_inversion(datetime.empty(0,1), ...
    'stride3_parameter_mean_2012_2015');
```

结果分别写入本目录的 `results/` 和 `figures/`，不会覆盖历史缓存。

已有全量结果时，可直接绘制模拟亮温减 SMOS 观测亮温的频率分布：

```matlab
cd('D:\海冰厚度联合反演\retrieval_matlab_ver\analysis\tb_spatial_averaging\baseline_parameter_mean_stride3')
plot_stride3_tb_residual_distribution();
```

差值定义为 `TB_sim - TB_SMOS`；正值表示模型偏暖。纵轴为每个区间的
样本百分比，不是概率密度。

继续拆解年份、冰型、OIB 覆盖度、飞行日及 `|差值|>=10 K` 长尾点：

```matlab
analyze_stride3_tb_residual_structure();
```

检查 `|差值|>=10 K` 长尾格的 SMOS `nPair`、亮温不确定度和 RFI：

```matlab
inspect_stride3_tail_smos_quality();
```
