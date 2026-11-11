# SMOS / SMAP 36 km EASE2-N maps

`run_ease2_sensor_maps_20160419` creates SMOS and SMAP 36 km
brightness-temperature maps for 2016-04-19 on the common EASE-Grid 2.0
North 12.5 km analysis grid (EPSG:6931, 1440 x 1440).

- SMOS: ICDC daily L3B `TB`, nearest-neighbour remap from EPSG:3413.
- SMAP 36 km: SPL1CTB.006, the same daily-composite rule.

The SMAP products are first averaged on their native EASE2-N grid and then
sampled at the target 12.5 km cell centres. This changes the analysis grid,
not the physical resolution. All panels use latitude >= 50 degrees and one
shared colour scale.

Expected input directories:

```text
data/smos/2016/SMOS_TB_north_20160419.nc
data/smap/SPL1CTB_006/2016/04/19/*.h5
```

Run in MATLAB:

```matlab
addpath('D:\海冰厚度联合反演\retrieval_matlab_ver\analysis\ease2_sensor_maps')
outputs = run_ease2_sensor_maps_20160419;       % SMOS + SMAP 36 km
```

The script writes 300 dpi PNG figures to `figures/` and statistics/MAT
data to `results/` under this directory.
