# 海冰厚度联合反演 MATLAB 主线

本仓库实现 OIB、SMOS/SMAP 与卫星干舷的空间配准、L 波段多层前向亮温、
冰厚/雪深联合反演，以及不同空间平均算子的对比试验。当前代码已经从
“2017 年论文纯复现”扩展为可配置的多源联合反演主线。

## 当前默认口径

- 空间网格：EASE2-N 12.5 km；历史 EASE 只能显式启用；
- OIB 窗口：3×3，中心抽样步长 3，相位 `[2 2]`，窗口互不重叠；
- 辐射传输：Burke 递推；`upperReflection` 仅作敏感性试验；
- 分层：20 层冰、1 层雪；
- 卤水体积分数：Cox–Weeks 原始值，进入 Vant 经验式前默认不截断到 70‰；
- 表面温度：默认固定 −30°C，可切换为逐点 OIB `KT19_surf`；
- 低温保护：默认进入卤水经验式的温度不低于 −30°C；KT19 冷端线性外推
  是独立实验分支，不会静默改变冻结基线。

以上参数只在 [`retrieval_config.m`](retrieval_config.m) 中定义。主流程不得
同时传入完整 `config` 和第二份 `operatorOptions`。

## 目录结构

| 路径 | 内容 |
| --- | --- |
| `+grid_resampling/` | EASE/EASE2 网格、经纬度落格和传感器重采样唯一实现 |
| `readdata/` | OIB、SMOS、SMAP 数据读取 |
| `plot/` | 配准、反演、残差和物理剖面绘图 |
| `analysis/` | 平均算子、KT19、SMOS–SMAP、分层等实验分支 |
| `tests/` | 单元测试、算子验证和全仓静态检查 |
| `过去的markdown/` | 已退出主线的历史说明，只用于追溯旧口径 |

网格函数必须显式使用 package 名，例如：

```matlab
grid = grid_resampling.ease_grid("EASE2");
smosEase = grid_resampling.regrid_smos_to_ease(smos, "EASE2");
[row, col] = grid_resampling.PSP2EASE(latitude, longitude, "EASE2");
```

工程根目录和 `readdata/` 不再保留这些函数的同名包装文件。

## 快速开始

需要 MATLAB R2023a 和 Mapping Toolbox。原始 OIB、SMOS、SMAP 数据不随仓库
提交，默认从仓库上一级的 `data/` 目录读取。

```matlab
cd('D:/海冰厚度联合反演/retrieval_matlab_ver')
addpath(genpath(pwd))

config = retrieval_config();
report = run_regression_suite();

main              % 构建同日空间配对
main_inversion    % 运行两条联合反演路线
```

完整配置支持受控覆盖：

```matlab
overrides = struct( ...
    'temperature', struct('source', "oib_kt19"), ...
    'model', struct('brineColdExtrapolationMode', ...
        "linear_below_minus30"));
config = retrieval_config(overrides);
```

未知字段会直接报错，避免拼写错误或脚本私自增加另一份固定参数。

## 当前验证状态

2026-11-11 的冻结验收结果：

- 34 个单元测试通过；
- 4 个空间平均/反演算子验证通过；
- 78 个 MATLAB 文件通过静态检查，0 条消息；
- 2012-03-22 和 2015-04-03 的真实 OIB–SMOS EASE2/KT19 调用链已跑通。

运行：

```matlab
report = run_regression_suite();
```

## 结果文件策略

仓库主线保存源码、测试、Markdown、Draw.io 和必要文档，不提交能够从脚本
重建的大批 `results/`、`figures/`、MAT、CSV、FIG、PNG、PDF、SVG 文件。
本地结果不会被清除，只是由 `.gitignore` 排除，避免每次实验把二进制结果
混进代码提交。

## 推荐 Git 推送流程

```powershell
git status --short --branch
git add -A
git diff --cached --stat
git commit -m "说明本次改动"
git fetch origin
git rebase origin/main
git push origin main
git ls-remote origin refs/heads/main
```

关键点是：`git commit -a` 不会包含新文件；推送前必须确认暂存区确实包含
新建源码，并先同步远端，避免 `non-fast-forward`。

## 后续任务入口

- [`analysis/多源亮温与卫星干舷后续任务交接.md`](analysis/多源亮温与卫星干舷后续任务交接.md)
- [`analysis/前向TB空间平均与冰雪分层后续任务交接.md`](analysis/前向TB空间平均与冰雪分层后续任务交接.md)
- [`analysis/tb_spatial_averaging/README.md`](analysis/tb_spatial_averaging/README.md)
- [`+grid_resampling/README.md`](+grid_resampling/README.md)
