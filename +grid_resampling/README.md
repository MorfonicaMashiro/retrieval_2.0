# 网格定义、坐标转换与重采样

本 MATLAB package 集中保存公共空间网格函数：

- `ease_grid`：建立原版 EASE 或 EASE2-N 12.5 km 网格定义；
- `PSP2EASE`：把经纬度点定位到目标网格行列；
- `regrid_smos_to_ease`：把规则 SMOS 极射赤面栅格最近邻映射到目标网格；
- `regrid_to_ease2_12p5`：统一分派 OIB、SMOS、SMAP 数据到 EASE2。

调用链：

```text
regrid_to_ease2_12p5
├─ OIB/SMAP 点记录 -> PSP2EASE -> 保留记录并附加行列号
└─ SMOS 栅格      -> regrid_smos_to_ease -> 最近邻目标栅格
                         └─ ease_grid
```

工程中只保留本目录这一份实现。所有调用均显式使用
`grid_resampling.函数名(...)`，避免 MATLAB 路径顺序选择到同名函数。
