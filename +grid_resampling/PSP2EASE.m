function [row, col, x, y] = PSP2EASE(latitude, longitude, version)
%PSP2EASE 将经纬度点定位到北极 EASE/EASE2 网格。
%   [ROW,COL,X,Y] = GRID_RESAMPLING.PSP2EASE(LAT,LON,VERSION)
%   先把经纬度投影为目标网格的米制 X/Y，再计算从 1 开始的矩阵行列号。
%
%   LAT/LON 可以是标量、向量或同形矩阵，例如 OIB 航迹点或传感器格心。
%   X/Y 保持输入形状，单位为 m；落在目标栅格外的 ROW/COL 设为 NaN，
%   但 X/Y 仍保留投影结果，便于诊断。
%
%   本函数只做“点落在哪个格子”的定位，不平均观测值，也不改变传感器
%   分辨率。SMOS 栅格最近邻重采样由 REGRID_SMOS_TO_EASE 负责。

    % 历史默认是原版 EASE；新流程显式传入 "EASE2"。
    if nargin < 3
        version = "EASE";
    end
    if ~isequal(size(latitude), size(longitude))
        error('PSP2EASE:SizeMismatch', ...
            'Latitude and longitude must have the same size.');
    end

    % 取得目标网格投影、边界和像元大小，然后把所有经纬度一次性投影。
    grid = grid_resampling.ease_grid(version);
    [x, y] = projfwd(grid.projection, ...
        double(latitude), double(longitude));

    % 行号从北向南增加：距北侧外边界 yMax 每跨过一个像元，行号加 1。
    % 列号从西向东增加：距西侧外边界 xMin 每跨过一个像元，列号加 1。
    % floor 而非 round 的原因是这里判断“落入哪个像元面”，不是寻找
    % 距离最近的格心。
    row = floor((grid.yMax - y) ./ grid.cellSize) + 1;
    col = floor((x - grid.xMin) ./ grid.cellSize) + 1;

    % 非有限投影结果或超出完整栅格边界的点不能形成合法矩阵索引。
    outside = ~isfinite(x) | ~isfinite(y) | ...
        row < 1 | row > grid.nRows | col < 1 | col > grid.nCols;
    row(outside) = NaN;
    col(outside) = NaN;
end
