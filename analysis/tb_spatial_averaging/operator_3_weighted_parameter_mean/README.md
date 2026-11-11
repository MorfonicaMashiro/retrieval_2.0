# 算子 3：到 SMOS 中心的反距离加权参数后模拟 TB

定义：

\[
TB_{Aw}=F\!\left(\frac{\sum_m w_mh_{i,m}}{\sum_m w_m},
                  \frac{\sum_m w_mh_{s,m}}{\sum_m w_m}\right).
\]

当前权重根据每个 OIB 点到对应 SMOS 中心格几何中心的距离计算：

\[
w_m=d_m^{-p},\qquad p=1.
\]

每个中心格分别计算其 3×3 邻域内的距离和权重，再归一化。若存在与中心重合
的点，则只对重合点等权，避免 `1/0`。`distancePower` 可改为 2 以加强近点权重。
这只是中心距离权重，不代表真实 SMOS 天线响应函数。

`collocate_day` 输出以下专属字段：

- `FBsnowWeighted`、`hsWeighted`、`FBiceWeighted`、`hiWeighted`
- `nOIBWeightedSupport`：参与该中心距离加权的 OIB 点数
- `TBOIBWeightedParameters`
- `TBWeightedMinusArithmetic`

正式加权参数反演入口：

```matlab
run_stride3_weighted_parameter_mean_inversion(datetime(2012,3,22), ...
    'stride3_distance_weighted_parameters_20120322')
```

该入口使用加权 `FBice/FBsnow` 求解，并以同一权重下的 `hi/hs` 作为验证值，
与本目录的加权参数前向算子保持一致。

入口：

```matlab
cd('D:\海冰厚度联合反演\retrieval_matlab_ver')
addpath('analysis/tb_spatial_averaging/operator_3_weighted_parameter_mean')
verify_weighted_parameter_mean_operator
run_weighted_parameter_mean_operator(datetime(2012,3,22), ...
    'stride3_distance_weighted_parameter_mean_20120322')
```

两个入口均固定使用 `centerStride=3`、`centerPhase=[2 2]`。步长 1 只存在于
已有历史结果中，后续不再生成。

本轮只改变参数的空间权重；表面温度仍固定为 −30 °C，使用 Burke 递推，
显式关闭 70‰ 卤水体积分数截断，且不计算逐 OIB 点 TB，避免混入算子 2 或
算子 4。目录中既有 `weighted_parameter_mean_*` 结果是旧的等面积子网格版本，
新的距离权重结果统一使用 `distance_weighted_*` 前缀，二者不会互相覆盖。

## 历史等面积版本

旧版等面积结果保留供对照，但不再作为本算子的默认权重。
