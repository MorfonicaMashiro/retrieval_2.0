# OIB KT19 表温 × EASE2 × 冷端线性外推

本分支不覆盖旧 EASE 冻结结果。新计算固定采用：

- EASE2-N 12.5 km（EPSG:6931）；
- 步长 3、相位 `[2 2]`、互不重叠 3×3 窗口；
- 先平均窗口内 OIB `hi/hs/KT19_surf`，再计算一次亮温；
- Burke 递推、20 层冰、1 层雪；
- Cox–Weeks 在 −30°C 以上使用原分段经验式；低于 −30°C 时，
  `F1(T)` 使用其在 −30°C 处的切线作连续线性外推；
- Vb 超过 70‰ 只标记，不截断。

全量运行：

```matlab
cd('D:\海冰厚度联合反演\retrieval_matlab_ver\analysis\tb_spatial_averaging\kt19_ease2_linear_extrapolation')
run_kt19_ease2_linear_extrapolation(datetime.empty(0,1), ...
    'kt19_ease2_linear_extrapolation_2012_2015');
```

在上述前向配对完成后，运行两条冰厚/雪深联合反演路线：

```matlab
run_kt19_ease2_inversion( ...
    'kt19_ease2_linear_extrapolation_2012_2015', ...
    'kt19_ease2_linear_extrapolation_2012_2015_inversion');
```

所有输出只写入本目录的 `results/` 与 `figures/`。代表日地图固定为
2012-03-22 和 2015-04-03，同时输出合并图和逐日图。

全量结果为 674 对：同样本固定 −30°C 的 RMSE 为 3.68 K，OIB KT19
输入的 RMSE 为 4.93 K。2015-03-24/25/26 因通过几何筛选的 OIB 点均
缺少有效 KT19 而无配对；2012-03-28 有 7,128 个有效 KT19 点，但与
当天有限 SMOS TB 的 EASE2 窗口交叠数为 0。逐日汇总表保存了
`tbOverlapCells`、`temperatureRejectedPoints` 和 `mixedTypeCells` 以区分原因。

联合反演中，TB+FBice 有解 483 对，冰厚/雪深 RMSE 分别为
0.496 m/0.169 m；TB+FBsnow 有解 470 对，冰厚/雪深 RMSE 分别为
0.811 m/0.125 m。在 `|ΔTB|≤1.5 K` 子集中，两条路线的冰厚 RMSE
分别为 0.306 m 和 0.276 m，雪深 RMSE 分别为 0.104 m 和 0.043 m。
