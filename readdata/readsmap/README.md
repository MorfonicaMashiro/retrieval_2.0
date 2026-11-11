# SMAP L1C 亮温读取

本目录提供两个正式入口：

- `read_smap_36km`：`SPL1CTB.006` 标准 36 km EASE-Grid 2.0；
- `read_smap_9km`：`SPL1CTB_E.004` 增强 9 km EASE-Grid 2.0。

两者共用 `read_smap_l1c_product`，只读取北极投影组。不会在读取阶段做日平均，
因为后续要使用逐格 `timeUtc` 与 OIB UTC 选择最近过境。

```matlab
matlabRoot = 'D:\海冰厚度联合反演\retrieval_matlab_ver';
addpath(matlabRoot, fullfile(matlabRoot, 'readdata', 'readsmap'));
config = retrieval_config();

data36 = read_smap_36km(file36, config.data.smap);
allObservations = data36.observations;
recommendedObservations = data36.recommended;

data9 = read_smap_9km(file9, config.data.smap);
```

`observations` 同时保留前视/后视以及下列信息：

- 原始与表面校正 H/V 亮温；
- 当前配置选择的 `tbH`、`tbV`；
- H/V 质量位和推荐质量掩膜；
- UTC、J2000 秒、入射角和扫描角；
- 网格行列、格心/观测质心坐标；
- H/V 样本数和水面比例；
- 产品、分辨率与源文件。

默认 `tbField="uncorrected"`，要求 H/V 的质量 bit 0 均为 0，亮温范围为
50--310 K。所有规则都位于 `retrieval_config().data.smap`。`recommended`
只是 `observations(recommendedMask,:)` 的便捷子表，原始记录不会丢失。

配置切换示例：

```matlab
overrides.data.brightnessTemperatureSource = "smap";
overrides.data.years = 2016:2019;
overrides.data.smap.activeProduct = "36km"; % 或 "9km"
config = retrieval_config(overrides);
```

当前主流程尚未把 SMAP 半轨合成和 OIB UTC 配准接入 `main.m`；后续应在读取后按
格点选择离 OIB 最近的过境，主分析使用 `|Δt| <= 3 h`，并报告 `<= 6 h`
敏感性。不要先把全部半轨做无条件日平均。
