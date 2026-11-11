function grid = ease_grid(version)
%EASE_GRID 构造北极 EASE/EASE2 12.5 km 网格的完整几何定义。
%   GRID = GRID_RESAMPLING.EASE_GRID(VERSION) 返回投影对象、行列数、
%   像元尺寸、左上外边界以及全部行列格心坐标。
%
%   VERSION="EASE"：历史球面 EASE-Grid North，1441×1441，
%   像元间距 12,533.76 m，名义 EPSG:3408。
%
%   VERSION="EASE2"：EASE-Grid 2.0 North，1440×1440，
%   像元间距 12,500 m，EPSG:6931。当前新分析默认使用该网格。
%
%   坐标约定：
%     x 向东增大；y 向北增大；矩阵行号从北向南增大；
%     xMin/yMax 表示整个栅格左上角的外边界，而不是首个像元中心。

    % 保留历史接口：未指定版本时返回原版 EASE。当前主程序会显式传入
    % "EASE2"，不会依赖这个历史默认值。
    if nargin < 1 || isempty(version)
        version = "EASE";
    end
    version = upper(string(version));

    % 两种网格不能只靠 EPSG 编号互换。原版 EASE 使用半径
    % 6,371,228 m 的球面；EASE2 使用标准 EPSG:6931 椭球投影。
    switch version
        case "EASE"
            n = 1441;
            cellSize = 12533.76;
            edge = 9030574.08;

            % MATLAB R2023a 中 projcrs(3408) 实际返回的是 EASE2 定义，
            % 因此这里手工建立原版球面 Lambert 方位等面积投影。
            projection = defaultm('eqaazim');
            projection.origin = [90 0 0];
            projection.geoid = [6371228 0];
            projection = defaultm(projection);
            epsg = 3408;

        case "EASE2"
            n = 1440;
            cellSize = 12500;
            edge = 9000000;
            projection = projcrs(6931);
            epsg = 6931;

        otherwise
            error('ease_grid:Version', ...
                '版本只能为 "EASE" 或 "EASE2"。');
    end

    % 先保存描述整个栅格的标量元数据。两种北极网格都是方阵，因此
    % nRows 与 nCols 相同。
    grid = struct('version', version, 'epsg', epsg, ...
        'nRows', n, 'nCols', n, ...
        'cellSize', cellSize, 'xMin', -edge, 'yMax', edge, ...
        'projection', projection);

    % 第一个格心位于外边界向内半个像元处。列号增加时 x 向东增加；
    % 行号增加时 y 向南减小。这两个一维数组用于定位和生成目标网格。
    grid.xCenters = grid.xMin + ((1:n) - 0.5) * cellSize;
    grid.yCenters = grid.yMax - ((1:n) - 0.5) * cellSize;
end
