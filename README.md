# 2017 论文数据配准（MATLAB）

本目录已完成 2017 年论文第 5.1 节的 **OIB—SMOS 同日配准、3×3 空间取样、四变量算术均值**。二分求解、理论解空间扫描和验证拟合的函数也已写好；辐射前向模型尚缺完整实现，当前不输出实际反演值。在 MATLAB R2023a 中进入本目录运行 `main`。需要 Mapping Toolbox。

## 论文所说的四个均值

2017 年论文 *Data Synergy between Altimetry and L-Band Passive Microwave Remote Sensing for the Retrieval of Sea Ice Parameters—A Theoretical Study of Methodology* 第 5.1 节说：对每次 OIB 飞行取 **同一天**的 SMOS TB；对每个 SMOS TB 中心格，把该格及周围八个 12.5 km EASE 格内的 OIB 测量作为样本，覆盖约 37.5 km × 37.5 km。排除混合冰型、开阔水和冰间水道后，对全部剩余 OIB **原始点**的以下四个量分别求算术均值：

| 论文记号 | 本地 OIB 字段或计算 | 单位 | 在论文验证中的作用 |
| --- | --- | --- | --- |
| `FBsnow` | `mean_fb`，雪面/总干舷 | m | 激光测高 + TB 路线的干舷输入 |
| `hs` | `snow_depth` | m | 雪深验证参考 |
| `FBice` | 逐点 `mean_fb - snow_depth`，冰面干舷代理 | m | 雷达测高 + TB 路线的干舷输入 |
| `hi` | `thickness` | m | 冰厚验证参考 |

若中心格 $j$ 的九格内有 $M_j$ 个通过筛选、四个量均有效的 OIB 点，则
$$
\overline q_j=\frac{1}{M_j}\sum_{p=1}^{M_j}q_p,
\qquad q\in\{FB_{snow},h_s,FB_{ice},h_i\}.
$$
同一批样本用于四个均值，因此 `mean(FBice) = mean(FBsnow) - mean(hs)`。这里按 **原始点数**加权；九个格子不会先各算一个均值再等权平均。邻近中心格的九格窗口会重叠，一个 OIB 点可能对应多个 SMOS TB，这正是按每个中心格取样的结果。

论文把 OIB 的 `FBice` 当作雷达干舷代理；它是雪面干舷减雪深得到的量，并非独立 CryoSat-2 观测。`hi` 和 `hs` 留作后续验证。目前脚本只存配对表，**不求反演值**。

## 坐标与配准顺序

```mermaid
flowchart LR
    A[同日 SMOS L3B<br/>原生极地立体网格] -->|最近邻重采样| C[12.5 km EASE 完整矩阵<br/>TB: 行, 列]
    B[OIB 航迹点<br/>经纬度与四变量] -->|投影及落格| D[OIB 所在 EASE 行, 列]
    C --> E[每个有效 TB 中心格]
    D -->|中心及周围 8 格<br/>筛选并对原始点求均值| E
    E --> F[同日配对表<br/>TB、四均值、nOIB]
```

1. `read_smos.m` 读取 ICDC SMOS L3B **逐日** NetCDF。原文件是 896×608 的北极极地立体网格，带逐像元经纬度；`read_oib.m` 读取 2012/2015 年 OIB 航次，并将数值缺测 `-99999` 转为 `NaN`。
2. `ease_grid.m` 定义完整北极原版 EASE 目标矩阵：1441×1441，步长 12,533.76 m，左上外边界 `(-9030574.08, 9030574.08)` m。像元中心为 `x(c)=xMin+(c-0.5)Δ`、`y(r)=yMax-(r-0.5)Δ`。
3. `regrid_smos_to_ease.m` 从每个目标中心逆向寻找最近的原 SMOS 像元，复制其 TB；原像元无 TB 或 `nPair<=0` 时保留缺测。只对 TB 做一次网格转换，**不对九个 SMOS TB 求均值**。
4. `PSP2EASE.m` 把每条 OIB 记录的经纬度投影到目标平面，按 `col=floor((x-xMin)/Δ)+1`、`row=floor((yMax-y)/Δ)+1` 定位行列；这是位置落格，不对 OIB 数值插值。
5. `collocate_2017_day.m` 先保留四变量均有效、`pcnt_ow=pcnt_thin_ice=pcnt_grey_ice=0`、`my_ice_flag` 为 FYI/MYI 的记录；再以每个 SMOS 中心格的 3×3 窗口汇总四个原始点总和与共同样本数。九格内同时含 FYI 与 MYI 的中心格被排除。第 5.1 节未要求负 `FBice` 阈值，因此此处不另加阈值。
6. `main.m` 顺序处理 2012 年 14 次与 2015 年 9 次 OIB 飞行，每次只取同日 SMOS。2012-03-28 有 7128 条 OIB 质量合格点，但对应九格区域内无有效 SMOS TB，因此该日无配对，详情看 `summary`。**不会**把 3–4 月亮温先求季节均值。

### PSP 与 EASE 是否有旋转

在当前代码采用的投影定义下，SMOS 源 PSP（EPSG:3413）的中央经线为 **−45°**，目标北极 EASE 的经度原点为 **0°**，所以同一经纬度在两套平面坐标中的极地方位相差 **45°**。但两者分别是极地立体投影和兰伯特方位等积投影，离极点的平面距离也不同，而且比例随纬度变化：例如 (80°N, 0°E) 在 PSP 中约为 (768, −768) km，在 EASE 中约为 (0, −1111) km。因此**不能把整张像素矩阵只旋转 45°**；代码通过目标格中心的经纬度在两个投影之间逐点转换，再取源 TB。

## 后续反演入口

2012 年 14 次飞行中，13 天有有效配对，共 3690 条；2015 年 9 天有配对，共 2852 条。论文的后续反演单位是**每一条配对**：先分别用 (TB, FBice) 和 (TB, FBsnow) 两条路线求解，再汇总多天的结果；不能先把 13 天 TB 或 OIB 值合成一个均值。

`main_inversion.m` 默认选择 2012 年并准备 `inputs` 和 `model`，可运行以检查记录数。`model` 含论文第 2.1 节的冬季温度、海水温盐和密度常数。`invert_2017_pair.m` 已实现论文图 5/6 的两条求解路线：雷达干舷以 4 cm 步长找区间再二分，雪面干舷以 1 cm 步长扫描所有变号区间再二分，并保留双解。论文没有规定亮温残差阈值，代码暂用 0.01 K。雪面干舷双解在第 5.2 节验证中默认取雪深较小的一支。

`solution_space_2017.m` 按第 2.2 节扫描 0–10 m 冰厚、0–1 m 雪深，依据静力平衡计算两种干舷并排除负冰面干舷；默认扫描间隔为 5 cm / 1 cm，**这是代码设置而非论文指定数值**。`plot/plot_solution_space_2017.m` 可绘制对应图 3/4 的 TB—干舷解空间。`fit_retrieval_2017.m` 计算图 11/12 使用的自由斜率最小二乘拟合、固定斜率为 1 的偏差和 RMSE。启用求解后，主程序也会按论文用模拟 OIB TB 与 SMOS TB 相差不超过 1.5 K 的样本另算一组验证统计。

`tb_forward_2017.m` 仍是式 (1) 的接口。论文 [32] 给出了多层结构、盐度剖面与混合像元公式，但没有在这两篇论文中完整列出可直接计算冰面辐射亮温的全部介电/辐射传输细节。`runSolver=false`，**不会输出实际反演值或理论解空间**；打开时若模型尚未补全，会明确报错。OIB 的 `hi` 和 `hs` 保留在 `inputs` 中供将来验证。

### 亮温前向模型的文献链

文献 [32] 的式 (5) 是海冰、水道、开阔水的**面积混合**，式 (6) 是把 3×3 EASE 区域内各 OIB 点的**已模拟海冰亮温**取均值。2017 联合反演只保留无水道、无开阔水的单一冰型样本，因而式 (5) 的面积权重可简化为 `αice=1`、其余为 0；但式 (6) 里的单点 `TBice(hi,hs,Tsfc,iceType,...)` 仍须计算。这里的 `TBice` 指整套雪—冰—海水层状系统在海冰覆盖区的辐射亮温，不是仅由海冰温度乘一个常数得到。

计算链为：候选冰厚、雪深和冰型 → 雪/冰各层温度与盐度 → 各层复介电常数 → 层间 Fresnel 反射、吸收和自发辐射 → 0–40°、两种极化的平均模拟 TB。2017 联合反演论文第 2.1 节把整条链压写成式 (1)。关键上游资料为：

1. [Maaß 等（2013）第 2.1 节](https://tc.copernicus.org/articles/7/1971/2013/)：明确采用 Burke 等（1979）的分层辐射传输模型，系统包含空气、雪、海冰与海水，说明温度剖面和材料介电常数的来源。
2. [Kaleschke 等（2010）第 2.2 节](https://tc.copernicus.org/articles/4/583/2010/)：给出冰的介电常数随盐水体积分数变化的公式，以及 1.4 GHz 的冰型系数表；其原始例子主要是**无雪单冰层**，不能直接代替当前多层雪冰模型。
3. [Zhou 等（2018）附录 A](https://tc.copernicus.org/articles/12/993/2018/)：同一作者对所用模型补充了 FYI/MYI 盐度剖面、热导率、雪和冰介电常数、盐水体积分数及顶层亮温式；其余光学/辐射传输细节继续引用前述模型及更早实验关系。

因此不必再找 OIB“密集度”补进当前反演。真正待实现的是 `TBice` 的逐层辐射计算，并对照论文图 2/3 与表 1 的模拟 TB 检查实现。公开的其他 Burke 模型实现可作编程参考，但不能未经核对就视为 Zhou 等（2017）的原作者代码。

**比较亮温所在层级：**SMOS 的 L1C 是大气层顶亮温，ICDC L3B 对 L1C 的水平/垂直极化与 0–40°入射角做日平均；[Kaleschke 等（2012）](https://doi.org/10.1029/2012GL050916)和[ICDC L3B 说明](https://www.cen.uni-hamburg.de/en/icdc/data/cryosphere/docs-cryo/documentation-smos-tbs.pdf)分别说明了这两步。Zhou 等（2018）附录 A1 写出 `TB_TOA`；在纯冰覆盖 `c=1` 时为 `TBice + (1-eice)*2.7 K + ΔTB_atm`。这里 `TBice` 是**海冰覆盖区的雪—冰—海水系统辐射分量**，不是冰表温度；海水仍是冰下的辐射边界。2017 年文献 [32] 的式 (5) 省写背景/大气项，并报告与 SMOS 比较的亮温偏差约 1.6 K。原论文固定的是水温、水盐度等输入参数，并未把 `TBice` 或 `TBwater` 固定为常数；附录 A1 也没有给出可直接计算逐日 `ΔTB_atm` 的完整方法。实现时必须明确采用哪一版亮温校正约定，不能把该偏差擅自解释成纯大气贡献。

### 密集度、密度与 OIB 字段

论文 [32] 的混合像元式 (5) 使用海冰、再冻结水道和开阔水的面积比例；这些比例来自独立海冰密集度产品和 MODIS 水道图，不是由 OIB 的雪深、干舷计算。2017 年联合反演论文第 2.1 节限制为无开阔水、无水道、单一冰型的冬季海冰；第 5.1 节排除混合区域。因此当前反演只需要 OIB 的冰型标记，以及 `FBice` / `FBsnow` 作为输入，`hi` / `hs` 用于验证。水、冰、雪密度固定为 1024/915/320 kg/m³，表温 −30°C、水温 −1.8°C、海水盐度 33 g/kg。

本地 OIB 表包含 `pcnt_ow`、`pcnt_thin_ice`、`pcnt_grey_ice`、`my_ice_flag`；当前配准用前三项及冰型标记筛除明显混合区域。它们不是论文 [32] 式 (5) 的独立海冰密集度与 MODIS 水道比例，不能把当前筛选视为完整水道掩膜。

若改做文献 [32] 的辐射模型验证，式 (6) 还会逐 OIB 点使用表面温度；本地原始表已有 `KT19_surf`，`read_oib.m` 读表时会保留它，但当前四变量均值配准不会聚合它。2017 联合反演论文的典型冬季实验直接采用固定 −30°C。

## 运行与输出

仓库包含当前的 `pairs_2017_same_day_v2.mat` 缓存，所以克隆后可直接运行 `main` 查看配准结果，运行 `main_inversion` 检查反演输入。若要从原始观测重新计算，把 `main.m` 的 `useCache` 改为 `false`，并将原始文件放在仓库同级的 `data/OIB` 和 `data/smos` 目录；原始观测数据不随仓库提交。

```matlab
cd('D:\海冰厚度联合反演\retrieval_matlab_ver')
main
head(pairs)
summary
main_inversion  % 查看逐条反演输入，目前不求物理解
```

- `processed/pairs_2017_same_day_v2.csv`：每条记录一个 SMOS 中心格，包含 `day, flight, row, col, TB, FBsnow, hs, FBice, hi, nOIB, isMYI`。
- `processed/pairs_2017_same_day_v2.mat`：同一配对表及各飞行日 `summary`；再次运行会直接读缓存。
- `figures/pairs_2017_coverage.png`：唯一的配准质检图，显示配对位置及九格内 OIB 样本数。

此前探索阶段生成的其他图片及缓存仍保留在 `figures/`、`processed/`，当前 `main` 不再使用。配准质检由 `plot/plot_2017_pairs.m` 绘制；补全前向模型后，`main_inversion` 才会调用 `plot/plot_solution_space_2017.m` 输出雷达和激光两种理论解空间图。

需要重算时在 `main.m` 将 `useCache = false`，运行后可改回 `true`。图像只绘制数千个有效配对点，避免为全部原始航迹反复生成大图。

## 复现边界

2017 年论文第 5.1 节称其 SMOS TB 为 **Level-4、12.5 km EASE 网格**，文献 [46] 指向 ICDC 的 *SMOS Daily Polar Gridded Brightness Temperatures*。本地文件也来自 ICDC 这一数据系列，但下载文件标为 **L3B v3、12.5 km 极地立体网格**。因此不是已证实换用了另一种 SMOS 观测产品；区别是论文所述的处理级别/最终网格与当前本地文件的发布标注/原始网格。当前代码把本地文件自己转换到 EASE，正是为了完成这一空间配准步骤。论文没有给出可唯一确认的 EASE 版本、重采样算法或其所用具体文件版本；本项目明确采用原版北极 EASE、EPSG:3413 源几何检查和最近邻 TB 重采样，**不能保证数值与作者完全一致**。原版 EASE 参数见 [NSIDC 指南](https://nsidc.org/data/user-resources/help-center/guide-ease-grids)，本地 SMOS 文件说明见 [ICDC L3B 数据页](https://www.cen.uni-hamburg.de/en/icdc/data/cryosphere/l3b-smos-tb.html)。

论文图 7、8 是理论解空间的输入误差敏感性图，并非这份 OIB—SMOS 配准表的可视化，当前阶段不绘制。用户提供的文献 [32] `Improving L-band radiation model and representation of small-scale variability to simulate brightness temperature of sea ice.pdf` 保留供后续讨论前向模型使用。
