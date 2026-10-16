function grid = ease_grid(version)
%EASE_GRID 建立 NSIDC 北极 12.5 km EASE 的完整目标栅格定义。
%   VERSION 为 "EASE" 或 "EASE2"。
%   xMin、yMax 是左上像元的外边界，xCenters、yCenters 是像元中心坐标。

    if nargin < 1 || isempty(version)
        version = "EASE";
    end
    switch upper(string(version))
        case "EASE"
            n = 1441; cellSize = 12533.76; edge = 9030574.08;
            % R2023a 的 projcrs(3408) 实际返回 EASE2，因此手动指定原版球面。
            projection = defaultm('eqaazim');
            projection.origin = [90 0 0];
            projection.geoid = [6371228 0];
            projection = defaultm(projection);
            epsg = 3408;
        case "EASE2"
            n = 1440; cellSize = 12500; edge = 9000000;
            projection = projcrs(6931);
            epsg = 6931;
        otherwise
            error('ease_grid:Version', '版本只能为 "EASE" 或 "EASE2"。');
    end

    grid = struct('version', upper(string(version)), 'epsg', epsg, ...
        'nRows', n, 'nCols', n, ...
        'cellSize', cellSize, 'xMin', -edge, 'yMax', edge, ...
        'projection', projection);
    grid.xCenters = grid.xMin + ((1:n) - 0.5) * cellSize;
    grid.yCenters = grid.yMax - ((1:n) - 0.5) * cellSize;
end
