# 算子 2 的第二轮温度试验：逐点 KT19 表面温度

本目录不修改固定 `surfaceTempC=-30 °C` 的冻结基线。它只在步长 3、
3×3 EASE 窗口、逐 OIB 点先算 TB 后平均的算子 2 上替换表面温度输入：

\[
TB_{KT19}=\frac{1}{M_T}\sum_{m\in\text{有效 KT19}}
F(h_{i,m},h_{s,m},T_{KT19,m},c_m).
\]

为避免把温度变化和样本缺测混在一起，结果同时计算完全相同 KT19 有效样本
在固定 −30 °C 下的 `TBOIBFixedMinus30SameKT19`。原算子 2 的全部样本固定
−30 °C 结果保存为 `TBOIBFixedMinus30All`。

`KT19_surf` 直接用于热传导、雪介电常数及各层辐射温度。由于 Cox–Weeks
卤水经验式的公开表值只到 −30 °C，本试验仅把进入该经验式的层温保护到
−30 °C；不把真实表面温度改回 −30 °C，也不在 −47 °C 处外推多项式。
对应配置为 `brineTemperatureFloorC=-30`，结果中会报告低于 −30/−35 °C 的
点数。Vb>70‰ 仍不截断，Burke 递推、20 层冰、1 层雪和其余参数保持不变。

单日测试：

```matlab
run_stride3_kt19_surface_temperature_experiment(datetime(2012,3,22), ...
    'kt19_surface_temperature_20120322')
```

全部 2012/2015 航次：

```matlab
run_stride3_kt19_surface_temperature_experiment(datetime.empty(0,1), ...
    'kt19_surface_temperature_2012_2015')
```

只输出 MAT、CSV 和彩色 PNG，不生成 SVG。

## 2012-03-22 单日验收结果

步长 3 得到 31 条配对，使用 28,968 个 OIB 点；这一天 `KT19_surf` 无缺测，
其中 28,259 点（97.55%）低于 −30 °C，14,618 点（50.47%）低于 −35 °C。

| 温度处理 | n | R² | Bias (K) | RMSE (K) |
|---|---:|---:|---:|---:|
| 固定 −30 °C，全部样本 | 31 | 0.7475 | −1.9267 | 2.1566 |
| 固定 −30 °C，同一 KT19 有效样本 | 31 | 0.7475 | −1.9267 | 2.1566 |
| 逐点 KT19，卤水经验式温度保护至 −30 °C | 31 | 0.7547 | −3.5460 | 3.6675 |

逐点 KT19 使模拟 TB 平均再降低约 1.62 K，相关趋势略有提高，但负偏差和
RMSE 明显增大。这只是单日结果；不能据此推断 2012/2015 全量表现。
