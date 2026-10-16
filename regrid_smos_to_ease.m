function result = regrid_smos_to_ease(smos, version)
%   REGRID_SMOS_TO_EASE 将 SMOS TB 最近邻重采样到北极 12.5 km EASE 网格。
%   使用逆向映射：每个目标 EASE 像元中心投到源 SMOS 极地网格，取最近源像元。
%   源 TB 缺测时结果仍为 NaN。VERSION 为 "EASE"（默认）或 "EASE2"，
%   目标像元名义分辨率均为 12.5 km。
%
%   最近邻是本项目明确选用的方法，不代表论文未公开的实际处理方法。
%   源网格几何会用文件自带的经纬度矩阵检查。

    if nargin < 2, version = "EASE"; end
    required = {'tb', 'latitude', 'longitude', 'nPair'};
    if ~all(isfield(smos, required))
        error('regrid_smos_to_ease:Fields', ...
            'Input must be a struct returned by read_smos.');
    end

    grid = ease_grid(version);
    sourceCRS = projcrs(3413);
    [sourceX, sourceY] = projfwd(sourceCRS, double(smos.latitude), double(smos.longitude));
    [nSourceRows, nSourceCols] = size(smos.tb);
    xAxis = sourceX(round(nSourceRows/2), :);
    yAxis = sourceY(:, round(nSourceCols/2));
    dx = median(diff(xAxis), 'omitnan');
    dy = median(diff(yAxis), 'omitnan');
    if ~isfinite(dx) || ~isfinite(dy) || abs(dx) < 1000 || abs(dy) < 1000
        error('regrid_smos_to_ease:GridGeometry', ...
            'Could not identify regular source polar-grid axes.');
    end
    residualX = abs(sourceX - xAxis);
    residualY = abs(sourceY - yAxis);
    residual = [residualX(isfinite(residualX)); residualY(isfinite(residualY))];
    if isempty(residual) || max(residual) > 0.05 * min(abs([dx, dy]))
        error('regrid_smos_to_ease:GridGeometry', ...
            'Source geolocation does not agree with a regular polar grid.');
    end

    [targetX, targetY] = meshgrid(grid.xCenters, grid.yCenters);
    [targetLat, targetLon] = projinv(grid.projection, targetX, targetY);
    [queryX, queryY] = projfwd(sourceCRS, targetLat, targetLon);
    clear targetX targetY targetLat targetLon

    colFloat = (queryX - xAxis(1)) ./ dx + 1;
    rowFloat = (queryY - yAxis(1)) ./ dy + 1;
    inside = isfinite(colFloat) & isfinite(rowFloat) & ...
        colFloat >= 0.5 & colFloat < nSourceCols + 0.5 & ...
        rowFloat >= 0.5 & rowFloat < nSourceRows + 0.5;
    sourceCol = round(colFloat(inside));
    sourceRow = round(rowFloat(inside));
    sourceIndices = sub2ind(size(smos.tb), sourceRow, sourceCol);

    tb = nan(grid.nRows, grid.nCols);
    nPair = nan(grid.nRows, grid.nCols);
    tb(inside) = smos.tb(sourceIndices);
    nPair(inside) = smos.nPair(sourceIndices);
    tb(~isfinite(nPair) | nPair <= 0) = NaN;

    result = struct('tb', tb, 'nPair', nPair, 'grid', grid, ...
        'sourceEPSG', 3413, 'method', "最近邻源像元");
end
