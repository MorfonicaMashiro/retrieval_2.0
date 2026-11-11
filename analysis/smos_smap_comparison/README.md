# 2015-04-01 / 2015-04-03 SMOS–SMAP 观测亮温对比

本目录是独立测试分支，不修改 2017 论文冻结基线，也不参与反演。

> **当前推荐口径已改为 SMAP 36 km 原生格共同支撑域。** SMOS 原版 EASE 12.53376 km 像元按真实多边形交叠面积汇总，并使用逐格 SMAP UTC 与 OIB 窗口时间比较。方法、时间诊断和新结果见 [`AREA_OVERLAP.md`](AREA_OVERLAP.md)。下文的“旧 EASE 格心落入 SMAP 格并作全天平均”结果只保留作历史对照。

## 比较口径

- 日期：2015-04-01、2015-04-03。
- 空间样本：沿用项目当前测试口径，只保留能够形成 OIB–SMOS 配对的步长 3 中心格，`centerPhase=[2 2]`。
- SMOS：ICDC L3B 北半球逐日 `TB`，先按现有程序重投影到原版 12.5 km EASE 网格。
- SMAP：NASA NSIDC DAAC `SPL1CTB` Version 6，36 km EASE-Grid 2.0；对所有通过主质量位的前视/后视观测分别求 H、V 日平均，再计算 `(H+V)/2`。
- 配准：把每个旧版 EASE 配对中心投影到 SMAP EASE2-N 36 km 网格，取所在 SMAP 网格。输出中心距离，超过半个 36 km 网格对角线的记录自动排除。
- 差值统一定义为 `SMAP - SMOS`。正值表示 SMAP 高于 SMOS。

## SMAP 产品等级与官方说明

本次下载的是 **Level-1C** 亮温产品，不是 L1B，也不是土壤水分等 L2/L3 地球物理反演产品：

- 产品短名：`SPL1CTB`
- 完整名称：SMAP L1C Radiometer Half-Orbit 36 km EASE-Grid Brightness Temperatures
- 版本：Version 6 (`SPL1CTB.006`)
- 内容：将已经定标、定位的 L1B 时序亮温重映射到 36 km EASE-Grid 2.0；每个文件对应半轨，包含约 1.41 GHz 的 H/V 极化和前视/后视亮温及质量标志。
- 官方产品页：https://nsidc.org/data/spl1ctb/versions/6
- 官方 User Guide（PDF）：https://nsidc.org/sites/default/files/documents/user-guide/spl1ctb-v006-userguide.pdf
- DOI：https://doi.org/10.5067/DV7IX2DQ681Y

官方说明指出，L1C 的核心处理是把 L1B 亮温映射到固定网格，并未进一步反演地球物理参数。因此它适合本次“观测亮温对观测亮温”的比较。

## 不能忽略的限制

这不是严格的传感器定标比较。当前 SMOS L3B `TB` 是论文采用的 H/V 与 0--40°综合量，而 SMAP 是约 40°固定入射角。`(H+V)/2` 只统一了极化平均，没有消除角度响应差异。因此本分支首先回答“两个现有观测产品在研究样本上的一致性如何”，不能仅凭差值判定任一传感器存在定标偏差。

官方给出的全局 EASE2 投影纬度边界约为 85.04°N；北极投影变量的实际覆盖及极区空洞仍应逐文件检查。本次被步长 3 选中的 68 个比较中心均获得了有效 SMAP 亮温，没有因覆盖范围造成缺测。

## 数据目录

把清单中的官方 HDF5 文件放入：

```text
data/smap/SPL1CTB_006/2015/04/01/
data/smap/SPL1CTB_006/2015/04/03/
```

所需文件及官方地址见 `smap_download_manifest.csv`。NASA 文件下载需要 Earthdata Login。

## MATLAB 命令

```matlab
addpath('D:\海冰厚度联合反演\retrieval_matlab_ver\analysis\smos_smap_comparison')
[pairs, metrics, summary] = run_smos_smap_two_day_comparison;
```

输出：

- `results/smos_smap_20150401_20150403_pairs.csv`
- `results/smos_smap_20150401_20150403_metrics.csv`
- `results/smos_smap_20150401_20150403_summary.csv`
- `results/smos_smap_20150401_20150403.mat`
- `figures/smos_smap_20150401_20150403_comparison.png`

## 本次结果

共得到 68 个唯一 SMAP 网格配对点（4 月 1 日 27 个，4 月 3 日 41 个），没有重复复用同一个 SMAP 网格。以 `(H+V)/2` 比较并定义差值为 `SMAP-SMOS`：

| 日期 | n | 相关系数 | R² | 平均差 / K | RMSE / K |
|---|---:|---:|---:|---:|---:|
| 2015-04-01 | 27 | 0.9231 | 0.8521 | -3.4741 | 3.8401 |
| 2015-04-03 | 41 | 0.8528 | 0.7273 | -2.9000 | 4.0023 |
| 两天合并 | 68 | 0.8719 | 0.7602 | -3.1279 | 3.9387 |

两天的平均差方向一致：在当前产品口径下，SMAP 的 H/V 平均亮温约比 SMOS 低 3.13 K。但由于入射角和产品合成方式尚未统一，这只能称为“产品间差值”，暂不能直接归因于 SMOS 或 SMAP 的系统定标偏差。
