function groups = build_stride3_point_groups(oib, pairs, grid, options)
%BUILD_STRIDE3_POINT_GROUPS 为互不重叠的 3x3 窗口保留逐 OIB 点观测。
%   本函数只组织数据，不计算 TB，也不使用 OIB 真值参与反演。

    if nargin < 4 || isempty(options), options = struct(); end
    defaults = struct('centerStride', 3, 'centerPhase', [2 2], ...
        'weightingMode', "uniform", 'subgridDivisions', 5, ...
        'excludeNegativeSnowFreeboard', true, ...
        'includeSurfaceTemperature', false, ...
        'surfaceTemperatureField', "KT19_surf", 'waterTempC', -1.8);
    names = fieldnames(defaults);
    for index = 1:numel(names)
        if ~isfield(options, names{index})
            options.(names{index}) = defaults.(names{index});
        end
    end
    options.weightingMode = string(options.weightingMode);
    if options.centerStride ~= 3 || ~isequal(options.centerPhase, [2 2])
        error('build_stride3_point_groups:Layout', ...
            '当前逐点分组只支持步长 3、相位 [2 2] 的互不重叠窗口。');
    end
    if ~ismember(options.weightingMode, ["uniform", "equal_area_subgrid"])
        error('build_stride3_point_groups:Weighting', ...
            'weightingMode 只能是 uniform 或 equal_area_subgrid。');
    end

    tbl = oib.table;
    required = {'lat','lon','mean_fb','snow_depth','thickness', ...
        'pcnt_ow','pcnt_thin_ice','pcnt_grey_ice','my_ice_flag'};
    if options.includeSurfaceTemperature
        required{end+1} = char(options.surfaceTemperatureField);
    end
    if ~all(ismember(required, tbl.Properties.VariableNames))
        error('build_stride3_point_groups:Fields', 'OIB 缺少逐点反演字段。');
    end
    fbSnow = double(tbl.mean_fb);
    hsTruth = double(tbl.snow_depth);
    hiTruth = double(tbl.thickness);
    fbIce = fbSnow - hsTruth;
    iceType = double(tbl.my_ice_flag);
    if options.includeSurfaceTemperature
        surfaceTempC = double(tbl.(char(options.surfaceTemperatureField)));
    else
        surfaceTempC = nan(height(tbl), 1);
    end
    valid = isfinite(tbl.lat) & isfinite(tbl.lon) & ...
        isfinite(fbSnow) & isfinite(hsTruth) & isfinite(hiTruth) & ...
        isfinite(fbIce) & tbl.pcnt_ow == 0 & ...
        tbl.pcnt_thin_ice == 0 & tbl.pcnt_grey_ice == 0 & ...
        ismember(iceType, [0 1]);
    if options.excludeNegativeSnowFreeboard
        valid = valid & fbSnow >= 0;
    end
    if options.includeSurfaceTemperature
        valid = valid & isfinite(surfaceTempC) & ...
            surfaceTempC < options.waterTempC;
    end
    [row, col, x, y] = grid_resampling.PSP2EASE( ...
        tbl.lat(valid), tbl.lon(valid), grid.version);
    located = isfinite(row) & isfinite(col) & isfinite(x) & isfinite(y);
    row = row(located); col = col(located);
    x = x(located); y = y(located);
    fbSnow = fbSnow(valid); fbSnow = fbSnow(located);
    hsTruth = hsTruth(valid); hsTruth = hsTruth(located);
    hiTruth = hiTruth(valid); hiTruth = hiTruth(located);
    fbIce = fbIce(valid); fbIce = fbIce(located);
    iceType = iceType(valid); iceType = iceType(located);
    surfaceTempC = surfaceTempC(valid); surfaceTempC = surfaceTempC(located);

    stride = options.centerStride;
    phase = options.centerPhase;
    centerRow = phase(1) + stride .* round((row - phase(1)) ./ stride);
    centerCol = phase(2) + stride .* round((col - phase(2)) ./ stride);
    inside = centerRow >= 1 & centerRow <= grid.nRows & ...
        centerCol >= 1 & centerCol <= grid.nCols;
    centerRow = centerRow(inside); centerCol = centerCol(inside);
    x = x(inside); y = y(inside);
    fbSnow = fbSnow(inside); hsTruth = hsTruth(inside);
    hiTruth = hiTruth(inside); fbIce = fbIce(inside);
    iceType = iceType(inside);
    surfaceTempC = surfaceTempC(inside);

    pointCenter = sub2ind([grid.nRows grid.nCols], centerRow, centerCol);
    pairCenter = sub2ind([grid.nRows grid.nCols], pairs.row, pairs.col);
    [belongs, pairIndex] = ismember(pointCenter, pairCenter);
    pairIndex = pairIndex(belongs);
    x = x(belongs); y = y(belongs);
    fbSnow = fbSnow(belongs); hsTruth = hsTruth(belongs);
    hiTruth = hiTruth(belongs); fbIce = fbIce(belongs);
    iceType = iceType(belongs);
    surfaceTempC = surfaceTempC(belongs);

    pointWeights = ones(size(fbSnow));
    if options.weightingMode == "equal_area_subgrid"
        fineCellSize = grid.cellSize / options.subgridDivisions;
        fineRows = grid.nRows * options.subgridDivisions;
        fineCols = grid.nCols * options.subgridDivisions;
        fineRow = floor((grid.yMax - y) ./ fineCellSize) + 1;
        fineCol = floor((x - grid.xMin) ./ fineCellSize) + 1;
        fineRow = min(max(fineRow, 1), fineRows);
        fineCol = min(max(fineCol, 1), fineCols);
        fineIndex = fineRow + (fineCol - 1) .* fineRows;
        combined = [double(pairIndex), double(fineIndex)];
        [~, ~, subgridGroup] = unique(combined, 'rows');
        countPerSubgrid = accumarray(subgridGroup, 1);
        pointWeights = 1 ./ countPerSubgrid(subgridGroup);
    end

    emptyGroup = struct('fbSnow', [], 'fbIce', [], 'hiTruth', [], ...
        'hsTruth', [], 'surfaceTempC', [], 'weights', [], 'isMYI', false, ...
        'nSamples', 0, 'nEffectiveCells', 0);
    groups = repmat(emptyGroup, height(pairs), 1);
    for currentPair = 1:height(pairs)
        selected = pairIndex == currentPair;
        if ~any(selected), continue; end
        currentType = iceType(selected);
        if any(currentType ~= currentType(1)) || ...
                logical(currentType(1)) ~= logical(pairs.isMYI(currentPair))
            error('build_stride3_point_groups:IceType', ...
                '逐点冰型与配对冰型不一致。');
        end
        currentWeights = pointWeights(selected);
        groups(currentPair).fbSnow = fbSnow(selected);
        groups(currentPair).fbIce = fbIce(selected);
        groups(currentPair).hiTruth = hiTruth(selected);
        groups(currentPair).hsTruth = hsTruth(selected);
        groups(currentPair).surfaceTempC = surfaceTempC(selected);
        groups(currentPair).weights = currentWeights;
        groups(currentPair).isMYI = logical(currentType(1));
        groups(currentPair).nSamples = nnz(selected);
        groups(currentPair).nEffectiveCells = sum(currentWeights);
    end
end
