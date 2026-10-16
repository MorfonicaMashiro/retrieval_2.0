function [row, col, x, y] = PSP2EASE(latitude, longitude, version)
%   PSP2EASE 将输入经纬度定位到北极 EASE 网格。
%   输入可为 OIB 航迹点或 SMOS 文件自带的经纬度矩阵。
%   X、Y 为目标投影的米坐标；
%   ROW、COL 为从 1 开始的 EASE 行列号。网格外点的行列号为 NaN。
%
%   VERSION 可取 "EASE"（原版，默认）或 "EASE2"；均为 12.5 km 网格。
%   本函数只转换位置，不重采样亮温；亮温栅格转换请调用 regrid_smos_to_ease。
%
%   OIB 直接用航迹经纬度；SMOS 可直接使用文件提供的经纬度。

    if nargin < 3, version = "EASE"; end
    if ~isequal(size(latitude), size(longitude))
        error('PSP2EASE:SizeMismatch', ...
            'Latitude and longitude must have the same size.');
    end

    grid = ease_grid(version);
    [x, y] = projfwd(grid.projection, double(latitude), double(longitude));
    row = floor((grid.yMax - y) ./ grid.cellSize) + 1;
    col = floor((x - grid.xMin) ./ grid.cellSize) + 1;
    outside = ~isfinite(x) | ~isfinite(y) | ...
        row < 1 | row > grid.nRows | col < 1 | col > grid.nCols;
    row(outside) = NaN;
    col(outside) = NaN;
end
