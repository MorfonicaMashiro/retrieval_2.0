# 算子 2：逐点干舷—逐点 TB 算术平均联合反演

正式路线在每个候选雪深下保留全部 OIB 干舷，由逐点干舷推导逐点冰厚，
逐点计算 TB 后再作足迹算术平均，并以该平均 TB 与 SMOS TB 求根。OIB
观测 `hi/hs` 只用于最终验证，不进入反演。

```matlab
cd('D:\海冰厚度联合反演\retrieval_matlab_ver')
addpath('analysis/tb_spatial_averaging/operator_2_sample_tb_mean')
verify_pointwise_freeboard_inversion
run_stride3_pointwise_freeboard_inversion(datetime(2012,3,22), ...
    'stride3_pointwise_uniform_20120322')
```

空日期数组会运行全部 2012/2015 航次。该计算会对每个候选雪深重新计算
全部局部 TB，显著慢于平均干舷路线；应先跑单日并核对结果。

## 旧的前向诊断

定义：

\[
TB_B=\frac{1}{M}\sum_{m=1}^{M}F(h_{i,m},h_{s,m}).
\]

以下旧脚本只复现“先逐样本前向、后平均”的验证口径。它与基线使用相同质量筛选、
冰型筛选、SMOS TB 和 3×3 EASE 邻域。现有 2012/2015 结果保存在
`results/`，历史结论见 `RESULTS_2012_2015.md`。它们复用了平均干舷反演根，
不能用于评价逐点反演的 `hi/hs` 精度。

入口：

```matlab
cd('D:\海冰厚度联合反演\retrieval_matlab_ver')
addpath('analysis/tb_spatial_averaging/operator_2_sample_tb_mean')
verify_sample_tb_mean_operator
run_sample_tb_mean_operator(datetime(2012,3,22), ...
    'tb_sample_mean_20120322')
```

新运行显式设置 `limitBrinePermittivity=false`。目录中已有结果生成于该设置
显式化之前，若用于最终论文比较，应按当前配置重新生成，不能与新结果混标。

步长 3、互不重叠窗口的 2012/2015 全量试验入口：

```matlab
cd('D:\海冰厚度联合反演\retrieval_matlab_ver\analysis\tb_spatial_averaging\operator_2_sample_tb_mean')
run_stride3_sample_tb_mean_experiment(datetime.empty(0,1), ...
    'stride3_sample_tb_mean_2012_2015');
```

该入口重新计算当前 Burke、无 70 截断设置下的逐 OIB 点 TB。反演输入没有
改变，因此它会先逐字段核对，再复用 `baseline_parameter_mean_stride3` 的求解根。
新结果分别写入 `results/stride3/` 与 `figures/stride3/`。
