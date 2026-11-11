# 算子 4：加权逐点干舷—逐点 TB 联合反演

定义：

\[
TB_{Bw}=\frac{\sum_m w_mF(h_{i,m},h_{s,m})}{\sum_m w_m}.
\]

该路线在每个候选雪深下，由逐点干舷推导逐点冰厚/雪深，逐点计算 TB 后，
按 2.5 km 等面积子网格权重汇总为足迹 TB，再与 SMOS TB 求根。权重用于
抑制 OIB 航迹采样密度，不代表真实 SMOS 天线响应。

```matlab
cd('D:\海冰厚度联合反演\retrieval_matlab_ver')
addpath('analysis/tb_spatial_averaging/operator_2_sample_tb_mean')
verify_pointwise_freeboard_inversion
addpath('analysis/tb_spatial_averaging/operator_4_weighted_sample_tb_mean')
run_stride3_weighted_pointwise_freeboard_inversion(datetime(2012,3,22), ...
    'stride3_pointwise_weighted_20120322')
```

结果独立写入 `results/pointwise_inversion/`。本路线不复用平均干舷反演根，
也不会覆盖算子 2 的等权结果。
