function result = regrid_smos_to_ease(smos, version)
%REGRID_SMOS_TO_EASE 将 SMOS 规则栅格最近邻映射到北极 EASE 网格。
%   RESULT = GRID_RESAMPLING.REGRID_SMOS_TO_EASE(SMOS,VERSION) 使用逆向
%   映射：对每个目标 EASE/EASE2 像元中心，查找投影坐标中最近的 SMOS
%   源像元中心，并复制其 TB、nPair 及可选质量字段。
%
%   SMOS 必须是 READ_SMOS 的输出。VERSION 可为 "EASE" 或 "EASE2"，
%   两种目标网格的名义间距约为 12.5 km。当前主线使用 "EASE2"。
%
%   重要限制：
%     * 这是格心最近邻，不是面积加权或双线性插值；
%     * 多个目标格可能复制同一个 SMOS 源 TB；
%     * 目标网格变细不会提高 SMOS 的真实物理分辨率；
%     * 最近源像元无效时结果为 NaN，不搜索第二近的有效像元。

    % 历史调用未传版本时保持原版 EASE；新分析应显式传 EASE2。
    if nargin < 2
        version = "EASE";
    end
    required = {'tb', 'latitude', 'longitude', 'nPair'};
    if ~all(isfield(smos, required))
        error('regrid_smos_to_ease:Fields', ...
            'Input must be a struct returned by read_smos.');
    end

    % 目标网格由 VERSION 决定。SMOS 文件本身位于北极极射赤面网格，
    % 因此先把文件自带的每个源像元经纬度转换到 EPSG:3413。
    grid = grid_resampling.ease_grid(version);
    sourceCRS = projcrs(3413);
    [sourceX, sourceY] = projfwd(sourceCRS, ...
        double(smos.latitude), double(smos.longitude));
    [nSourceRows, nSourceCols] = size(smos.tb);

    % SMOS 源数据应为与 EPSG:3413 坐标轴平行的规则矩形网格：
    % 同一列的 x 近似不变，同一行的 y 近似不变。取中间一行的 x 和
    % 中间一列的 y，恢复一维源网格轴；中部通常比边缘更完整。
    middleRow = round(nSourceRows / 2);
    middleCol = round(nSourceCols / 2);
    xAxis = sourceX(middleRow, :);       % 每个源网格列中心的 x
    yAxis = sourceY(:, middleCol);       % 每个源网格行中心的 y

    % 相邻格心坐标差即源网格间距。用中位数抑制个别定位异常和投影
    % 浮点误差。dx/dy 的负号仅表示文件中行列的排列方向。
    dx = median(diff(xAxis), 'omitnan');
    dy = median(diff(yAxis), 'omitnan');
    if ~isfinite(dx) || ~isfinite(dy) || ...
            abs(dx) < 1000 || abs(dy) < 1000
        error('regrid_smos_to_ease:GridGeometry', ...
            'Could not identify regular source polar-grid axes.');
    end

    % 检验规则矩形网格假设：理论上 sourceX 只由列号决定，sourceY 只
    % 由行号决定。若二维定位与一维坐标轴的偏差超过网格间距的 5%，
    % 就拒绝继续，避免错误使用后面的快速浮点行列号公式。
    residualX = abs(sourceX - xAxis);
    residualY = abs(sourceY - yAxis);
    residual = [residualX(isfinite(residualX)); ...
        residualY(isfinite(residualY))];
    if isempty(residual) || ...
            max(residual) > 0.05 * min(abs([dx, dy]))
        error('regrid_smos_to_ease:GridGeometry', ...
            'Source geolocation does not agree with a regular polar grid.');
    end

    % 生成全部目标像元中心，并依次转换为：目标投影坐标 -> 经纬度 ->
    % SMOS 源投影坐标。queryX/queryY 表示每个目标格心在源坐标系中的
    % 位置；这里采用逆向映射可避免源像元散投影后产生空洞或重复竞争。
    [targetX, targetY] = meshgrid(grid.xCenters, grid.yCenters);
    [targetLat, targetLon] = projinv( ...
        grid.projection, targetX, targetY);
    [queryX, queryY] = projfwd(sourceCRS, targetLat, targetLon);
    clear targetX targetY targetLat targetLon

    % 把目标格心的源投影坐标换算为 SMOS 浮点列号/行号。例如
    % colFloat=123.27、rowFloat=86.72 表示位于相邻源格心之间。
    colFloat = (queryX - xAxis(1)) ./ dx + 1;
    rowFloat = (queryY - yAxis(1)) ./ dy + 1;

    % 0.5 与 N+0.5 是首末源像元中心向外各延伸半个像元得到的栅格
    % 外边界。只有位于边界内的目标格心才允许建立源索引。
    inside = isfinite(colFloat) & isfinite(rowFloat) & ...
        colFloat >= 0.5 & colFloat < nSourceCols + 0.5 & ...
        rowFloat >= 0.5 & rowFloat < nSourceRows + 0.5;

    % 四舍五入浮点行列号，等价于选择规则投影坐标中最近的源格心。
    % 这里不是遍历所有源格计算球面距离，但对规则轴对齐网格结果相同。
    sourceCol = round(colFloat(inside));
    sourceRow = round(rowFloat(inside));
    sourceIndices = sub2ind( ...
        size(smos.tb), sourceRow, sourceCol);

    % 用同一最近邻索引复制 TB 与有效配对数。最近像元的 nPair 无效或
    % 不大于 0 时，目标 TB 直接设为 NaN，不继续寻找其他有效像元。
    tb = sample_source_raster(smos.tb, inside, sourceIndices, grid);
    nPair = sample_source_raster( ...
        smos.nPair, inside, sourceIndices, grid);
    tb(~isfinite(nPair) | nPair <= 0) = NaN;

    % 保存重采样结果和方法元数据，明确区分“目标网格间距”与“传感器
    % 物理分辨率”。sourceGridSpacingKm 是源文件规则网格的格心间距。
    result = struct('tb', tb, 'nPair', nPair, 'grid', grid, ...
        'sourceEPSG', 3413, ...
        'sourceGridSpacingKm', mean(abs([dx, dy])) / 1000, ...
        'method', "最近邻源像元");

    % 质量字段采用与 TB 完全相同的最近邻索引，确保每个目标格的 TB、
    % 不确定度和 RFI 信息仍来自同一个源像元。
    if isfield(smos, 'tbUncertainty')
        result.tbUncertainty = sample_source_raster( ...
            smos.tbUncertainty, inside, sourceIndices, grid);
    end
    if isfield(smos, 'rfiRatio')
        result.rfiRatio = sample_source_raster( ...
            smos.rfiRatio, inside, sourceIndices, grid);
    end
    if isfield(smos, 'source')
        result.source = smos.source;
    end
    if isfield(smos, 'day')
        result.day = smos.day;
    end
    result.resolutionNote = ...
        "12.5 km is target spacing; nearest-neighbor mapping does not add information.";
end

function target = sample_source_raster( ...
        source, inside, sourceIndices, grid)
%SAMPLE_SOURCE_RASTER 按已计算的最近邻索引复制一幅源栅格。
%   目标网格先填 NaN，只在落入源栅格边界的目标位置写入源值。因此源
%   覆盖范围以外的格子自然保持 NaN，且不同字段共享完全相同的几何映射。

    target = nan(grid.nRows, grid.nCols);
    target(inside) = double(source(sourceIndices));
end
