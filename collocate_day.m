function [pairs, stats] = collocate_day( ...
        oib, smosEase, day, configOrModel, forwardFcn, operatorOptions)
%COLLOCATE_DAY 汇集同日 OIB 与网格化被动微波亮温。
%   [PAIRS,STATS] = COLLOCATE_DAY(OIB,SMOSEASE,DAY,CONFIG) 是推荐接口。
%   CONFIG 同时提供前向模型、配准算子和 fixed/oib_kt19 表温选择。
%   为兼容已有实验，仍接受旧式 MODEL、FORWARDFCN、OPERATOROPTIONS。
%   每个有效 SMOS 中心格使用周围 3×3 EASE 格内的 OIB 原始记录。
%   对同一批有效记录分别求 FBsnow、hs、FBice、hi 的算术均值。
%   FBice = FBsnow - hs，是 OIB 几何冰面干舷代理。
%   本地质量字段中 0 表示无开阔水/薄冰/灰冰，my_ice_flag 为
%   0=FYI、1=MYI；九格内同时出现两种冰型的中心格会被排除。
%   可选传入 MODEL 后，还会对同一批 OIB 点逐点调用 FORWARDFCN，再在
%   3×3 邻域内求亮温算术平均。FORWARDFCN 默认是 @tb_forward。
%   OPERATOROPTIONS.parameterWeighting="equal_area_subgrid" 时，先把
%   OIB 点投影到更细的等面积子网格；同一子网格内的全部点合计权重为 1，
%   再对 3×3 邻域内已占用子网格等权，得到加权 hi/hs 后计算一次 TB。
%   parameterWeighting="inverse_distance_to_center" 时，对每个 SMOS 中心格
%   单独计算 OIB 点到该格中心的距离，使用 w=1/d^distancePower 加权。
%   若点与中心重合，则只对所有重合点等权，避免距离奇点。
%   OPERATOROPTIONS.subgridDivisions 设置每个 12.5 km EASE 格沿单轴细分数，
%   默认 5，即 2.5 km 子网格。computeSampleMeanTB=false 可只计算加权参数
%   算子，跳过逐 OIB 点前向计算。
%   OPERATOROPTIONS.windowSize 设置正方形邻域边长，论文基线为 3。
%   OPERATOROPTIONS.centerStride 设置 SMOS 中心格步长，当前默认 3。
%   centerStride=3、centerPhase=[2 2] 时只保留 (2,2)、(2,5)... 等中心，
%   对应从 EASE 网格起点划分的互不重叠 3×3 OIB 窗口。
%   原三参数调用不计算前向亮温，输出字段与历史版本保持一致。

    if nargin < 4, configOrModel = []; end
    isConfigInterface = isstruct(configOrModel) && ...
        all(isfield(configOrModel, {'model','collocation','temperature'}));
    if isConfigInterface && nargin > 4
        error('collocate_day:DuplicateConfiguration', ...
            ['完整 CONFIG 已包含前向模型和配准选项；不要再传入 ' ...
            'FORWARDFCN 或 OPERATOROPTIONS。']);
    end
    temperatureConfig = struct('source', "fixed", ...
        'fixedC', -30, 'oibField', "KT19_surf", ...
        'oibBrineFloorC', -30);
    if isConfigInterface
        config = configOrModel;
        model = config.model;
        forwardFcn = @tb_forward;
        operatorOptions = config.collocation;
        temperatureConfig = config.temperature;
        hasForwardModel = true;
    else
        model = configOrModel;
        hasForwardModel = ~isempty(model);
        if nargin < 5 || isempty(forwardFcn), forwardFcn = @tb_forward; end
        if nargin < 6 || isempty(operatorOptions), operatorOptions = struct(); end
        if hasForwardModel && isfield(model, 'surfaceTempC')
            temperatureConfig.fixedC = model.surfaceTempC;
        end
    end
    temperatureSource = string(temperatureConfig.source);
    usesOibTemperature = hasForwardModel && temperatureSource == "oib_kt19";
    operatorOptions = validate_operator_options( ...
        operatorOptions, hasForwardModel);
    computeSampleMeanTB = hasForwardModel && ...
        operatorOptions.computeSampleMeanTB;
    computeWeightedParameters = ...
        operatorOptions.parameterWeighting ~= "none";

    tbl = oib.table;
    required = {'pcnt_ow', 'pcnt_thin_ice', 'pcnt_grey_ice', 'my_ice_flag'};
    if ~all(ismember(required, tbl.Properties.VariableNames))
        error('collocate_day:Fields', 'OIB 缺少冰型或混合像元字段。');
    end
    grid = smosEase.grid;
    fbSnow = double(tbl.mean_fb);
    hs = double(tbl.snow_depth);
    fbIce = fbSnow - hs;
    hi = double(tbl.thickness);
    iceType = double(tbl.my_ice_flag);
    surfaceTemperature = repmat(double(temperatureConfig.fixedC), ...
        height(tbl), 1);
    if usesOibTemperature
        temperatureField = char(temperatureConfig.oibField);
        if ~ismember(temperatureField, tbl.Properties.VariableNames)
            error('collocate_day:TemperatureField', ...
                'OIB 缺少配置指定的表温字段 %s。', temperatureField);
        end
        surfaceTemperature = double(tbl.(temperatureField));
    end

    % 论文要求排除混合冰型、开阔水和冰间水道；这里使用产品质量字段落实。
    % 论文第 5.1 节没有要求排除负冰面干舷，因此不额外设置该阈值。
    geometryValid = isfinite(tbl.lat) & isfinite(tbl.lon) & ...
        isfinite(fbSnow) & isfinite(hs) & isfinite(hi) & ...
        isfinite(fbIce) & ...
        tbl.pcnt_ow == 0 & tbl.pcnt_thin_ice == 0 & ...
        tbl.pcnt_grey_ice == 0 & ismember(iceType, [0 1]);
    temperatureValid = true(height(tbl), 1);
    if usesOibTemperature
        temperatureValid = isfinite(surfaceTemperature) & ...
            surfaceTemperature < model.waterTempC;
    end
    valid = geometryValid & temperatureValid;
    [row, col, x, y] = grid_resampling.PSP2EASE( ...
        tbl.lat(valid), tbl.lon(valid), grid.version);
    located = isfinite(row) & isfinite(col);
    row = row(located);
    col = col(located);
    x = x(located);
    y = y(located);
    fbSnow = fbSnow(valid); fbSnow = fbSnow(located);
    hs = hs(valid); hs = hs(located);
    fbIce = fbIce(valid); fbIce = fbIce(located);
    hi = hi(valid); hi = hi(located);
    iceType = iceType(valid); iceType = iceType(located);
    surfaceTemperature = surfaceTemperature(valid);
    surfaceTemperature = surfaceTemperature(located);

    shape = [grid.nRows, grid.nCols];
    nCells = prod(shape);
    index = sub2ind(shape, row, col);
    count = accumarray(index, 1, [nCells 1], @sum, 0);
    fyiCount = accumarray(index, double(iceType == 0), [nCells 1], @sum, 0);
    myiCount = accumarray(index, double(iceType == 1), [nCells 1], @sum, 0);
    fields = {fbSnow, hs, fbIce, hi};
    surfaceTemperatureSum = [];
    if usesOibTemperature
        cellTemperatureSum = accumarray( ...
            index, surfaceTemperature, [nCells 1], @sum, 0);
        surfaceTemperatureSum = conv2(reshape( ...
            cellTemperatureSum, shape), ones(operatorOptions.windowSize), ...
            'same');
    end
    pointForwardSeconds = NaN;
    if computeSampleMeanTB
        forwardTimer = tic;
        pointTB = nan(size(hi));
        isFYI = iceType == 0;
        isMYI = iceType == 1;
        fyiModel = model;
        myiModel = model;
        if usesOibTemperature
            fyiModel.surfaceTempC = surfaceTemperature(isFYI);
            myiModel.surfaceTempC = surfaceTemperature(isMYI);
        end
        pointTB(isFYI) = forwardFcn( ...
            hi(isFYI), hs(isFYI), false, fyiModel);
        pointTB(isMYI) = forwardFcn( ...
            hi(isMYI), hs(isMYI), true, myiModel);
        if any(~isfinite(pointTB))
            error('collocate_day:ForwardTB', ...
                '逐 OIB 点前向亮温包含非有限值。');
        end
        fields{end + 1} = pointTB;
        pointForwardSeconds = toc(forwardTimer);
    end
    sums = cell(1, numel(fields));
    kernel = ones(operatorOptions.windowSize);
    for k = 1:numel(fields)
        cellSum = accumarray(index, fields{k}, [nCells 1], @sum, 0);
        sums{k} = conv2(reshape(cellSum, shape), kernel, 'same');
    end
    nOIB = conv2(reshape(count, shape), kernel, 'same');
    nFYI = conv2(reshape(fyiCount, shape), kernel, 'same');
    nMYI = conv2(reshape(myiCount, shape), kernel, 'same');
    weightedSums = {};
    weightedDenominator = [];
    weightedSupport = [];
    if computeWeightedParameters
        if operatorOptions.parameterWeighting == "equal_area_subgrid"
            pointWeights = equal_area_subgrid_weights( ...
                x, y, grid, operatorOptions.subgridDivisions);
            cellWeight = accumarray( ...
                index, pointWeights, [nCells 1], @sum, 0);
            weightedDenominator = conv2( ...
                reshape(cellWeight, shape), kernel, 'same');
            weightedSupport = weightedDenominator;
            weightedFields = fields(1:4);
            if usesOibTemperature
                weightedFields{end + 1} = surfaceTemperature;
            end
            weightedSums = cell(1, numel(weightedFields));
            for k = 1:numel(weightedFields)
                cellWeightedSum = accumarray(index, ...
                    pointWeights .* weightedFields{k}, [nCells 1], @sum, 0);
                weightedSums{k} = conv2( ...
                    reshape(cellWeightedSum, shape), kernel, 'same');
            end
        else
            weightedFields = fields(1:4);
            if usesOibTemperature
                weightedFields{end + 1} = surfaceTemperature;
            end
            [weightedSums, weightedDenominator, weightedSupport] = ...
                inverse_distance_center_accumulators( ...
                row, col, x, y, weightedFields, grid, ...
                operatorOptions.distancePower, operatorOptions.windowSize);
        end
    end
    good = nOIB > 0 & isfinite(smosEase.tb) & ...
        ((nFYI > 0 & nMYI == 0) | (nMYI > 0 & nFYI == 0));
    if operatorOptions.centerStride > 1
        selectedRows = mod((1:shape(1)).' - ...
            operatorOptions.centerPhase(1), ...
            operatorOptions.centerStride) == 0;
        selectedCols = mod((1:shape(2)) - ...
            operatorOptions.centerPhase(2), ...
            operatorOptions.centerStride) == 0;
        good = good & selectedRows & selectedCols;
    end
    centers = find(good);
    [centerRow, centerCol] = ind2sub(shape, centers);
    sampleCount = nOIB(centers);
    [~, flightName] = fileparts(char(oib.source));
    n = numel(centers);
    pairs = table(repmat(day, n, 1), repmat(string(flightName), n, 1), ...
        centerRow, centerCol, smosEase.tb(centers), ...
        sums{1}(centers) ./ sampleCount, ...
        sums{2}(centers) ./ sampleCount, ...
        sums{3}(centers) ./ sampleCount, ...
        sums{4}(centers) ./ sampleCount, sampleCount, nMYI(centers) > 0, ...
        'VariableNames', {'day','flight','row','col','TB', ...
        'FBsnow','hs','FBice','hi','nOIB','isMYI'});
    pairForwardTimer = tic;
    forwardPairEvaluations = 0;
    if hasForwardModel
        centerIsMYI = nMYI(centers) > 0;
        centerIsFYI = ~centerIsMYI;
        pairModel = model;
        if usesOibTemperature
            meanSurfaceTemperature = ...
                surfaceTemperatureSum(centers) ./ sampleCount;
            pairModel.surfaceTempC = meanSurfaceTemperature;
            pairs = addvars(pairs, meanSurfaceTemperature, ...
                'NewVariableNames', 'surfaceTempC');
        end
        tbOIBMeanParameters = evaluate_pair_forward( ...
            pairs.hi, pairs.hs, centerIsMYI, centerIsFYI, ...
            pairModel, forwardFcn);
        forwardPairEvaluations = forwardPairEvaluations + n;
        pairs = addvars(pairs, tbOIBMeanParameters, ...
            'NewVariableNames', 'TBOIBMeanParameters');
    end
    if computeSampleMeanTB
        tbOIBMeanSamples = sums{5}(centers) ./ sampleCount;
        tbSampleMinusParameter = ...
            tbOIBMeanSamples - pairs.TBOIBMeanParameters;
        pairs = addvars(pairs, tbOIBMeanSamples, ...
            tbSampleMinusParameter, 'NewVariableNames', ...
            {'TBOIBMeanSamples','TBSampleMinusParameter'});
    end
    if computeWeightedParameters
        weightSum = weightedDenominator(centers);
        weightedValues = cellfun( ...
            @(value) value(centers) ./ weightSum, weightedSums, ...
            'UniformOutput', false);
        fbSnowWeighted = weightedValues{1};
        hsWeighted = weightedValues{2};
        fbIceWeighted = weightedValues{3};
        hiWeighted = weightedValues{4};
        weightedModel = model;
        if usesOibTemperature
            surfaceTempCWeighted = weightedValues{5};
            weightedModel.surfaceTempC = surfaceTempCWeighted;
        end
        nOIBWeightedSupport = round(weightedSupport(centers));
        if operatorOptions.parameterWeighting == "equal_area_subgrid"
            nOIBSubgrid = nOIBWeightedSupport;
        else
            nOIBSubgrid = nan(size(nOIBWeightedSupport));
        end
        tbOIBWeightedParameters = evaluate_pair_forward( ...
            hiWeighted, hsWeighted, centerIsMYI, centerIsFYI, ...
            weightedModel, forwardFcn);
        forwardPairEvaluations = forwardPairEvaluations + n;
        tbWeightedMinusArithmetic = ...
            tbOIBWeightedParameters - pairs.TBOIBMeanParameters;
        pairs = addvars(pairs, fbSnowWeighted, hsWeighted, ...
            fbIceWeighted, hiWeighted, nOIBWeightedSupport, nOIBSubgrid, ...
            tbOIBWeightedParameters, tbWeightedMinusArithmetic, ...
            'NewVariableNames', {'FBsnowWeighted','hsWeighted', ...
            'FBiceWeighted','hiWeighted','nOIBWeightedSupport', ...
            'nOIBSubgrid', ...
            'TBOIBWeightedParameters', ...
            'TBWeightedMinusArithmetic'});
        if usesOibTemperature
            pairs = addvars(pairs, surfaceTempCWeighted, ...
                'NewVariableNames', 'surfaceTempCWeighted');
        end
    end
    pairForwardSeconds = toc(pairForwardTimer);
    if ~hasForwardModel
        pairForwardSeconds = NaN;
    end
    forwardSeconds = sum([pointForwardSeconds, pairForwardSeconds], ...
        'omitnan');
    if ~hasForwardModel
        forwardSeconds = NaN;
    end
    stats = struct('rawPoints', height(tbl), ...
        'geometryQualityPoints', nnz(geometryValid), ...
        'temperatureRejectedPoints', nnz(geometryValid & ~temperatureValid), ...
        'qualityPoints', numel(index), ...
        'tbOverlapCells', nnz(nOIB > 0 & isfinite(smosEase.tb)), ...
        'matchedCells', n, 'mixedTypeCells', ...
        nnz(nOIB > 0 & nFYI > 0 & nMYI > 0 & isfinite(smosEase.tb)), ...
        'forwardPoints', double(computeSampleMeanTB) * numel(index), ...
        'forwardPairs', forwardPairEvaluations, ...
        'pointForwardSeconds', pointForwardSeconds, ...
        'pairForwardSeconds', pairForwardSeconds, ...
        'forwardSeconds', forwardSeconds, ...
        'parameterWeighting', operatorOptions.parameterWeighting, ...
        'subgridDivisions', operatorOptions.subgridDivisions, ...
        'distancePower', operatorOptions.distancePower, ...
        'centerStride', operatorOptions.centerStride, ...
        'centerPhase', operatorOptions.centerPhase, ...
        'surfaceTemperatureSource', temperatureSource);
end

function options = validate_operator_options(options, hasForwardModel)
    persistent baselineDefaults
    if isempty(baselineDefaults)
        config = retrieval_config();
        baselineDefaults = config.collocation;
    end
    defaults = baselineDefaults;
    % 通用接口保持历史行为：传入 MODEL 时默认计算逐点 TB；若只覆盖
    % centerStride 而未给 phase，则按步长自动选中心相位。
    defaults.computeSampleMeanTB = hasForwardModel;
    defaults.centerPhase = [];
    names = fieldnames(options);
    unknown = setdiff(names, fieldnames(defaults));
    if ~isempty(unknown)
        error('collocate_day:OperatorOptions', ...
            '未知算子选项：%s。', strjoin(unknown, ', '));
    end
    for k = 1:numel(names)
        defaults.(names{k}) = options.(names{k});
    end
    options = defaults;
    if ~(isscalar(options.computeSampleMeanTB) && ...
            (islogical(options.computeSampleMeanTB) || ...
            ismember(options.computeSampleMeanTB, [0 1])))
        error('collocate_day:OperatorOptions', ...
            'computeSampleMeanTB 必须是逻辑标量。');
    end
    options.computeSampleMeanTB = logical(options.computeSampleMeanTB);
    if ~(isscalar(options.windowSize) && isnumeric(options.windowSize) && ...
            isfinite(options.windowSize) && options.windowSize >= 1 && ...
            options.windowSize == fix(options.windowSize) && ...
            mod(options.windowSize, 2) == 1)
        error('collocate_day:OperatorOptions', ...
            'windowSize 必须是正奇数。');
    end
    options.parameterWeighting = string(options.parameterWeighting);
    if ~isscalar(options.parameterWeighting) || ...
            ~ismember(options.parameterWeighting, ...
            ["none", "equal_area_subgrid", ...
            "inverse_distance_to_center"])
        error('collocate_day:OperatorOptions', ...
            ['parameterWeighting 只能是 none、equal_area_subgrid 或 ', ...
            'inverse_distance_to_center。']);
    end
    if ~(isscalar(options.subgridDivisions) && ...
            isfinite(options.subgridDivisions) && ...
            options.subgridDivisions >= 1 && ...
            options.subgridDivisions == fix(options.subgridDivisions))
        error('collocate_day:OperatorOptions', ...
            'subgridDivisions 必须是正整数。');
    end
    if ~(isscalar(options.distancePower) && ...
            isnumeric(options.distancePower) && ...
            isfinite(options.distancePower) && options.distancePower > 0)
        error('collocate_day:OperatorOptions', ...
            'distancePower 必须是有限正数。');
    end
    if ~(isscalar(options.centerStride) && ...
            isnumeric(options.centerStride) && ...
            isfinite(options.centerStride) && ...
            options.centerStride >= 1 && ...
            options.centerStride == fix(options.centerStride))
        error('collocate_day:OperatorOptions', ...
            'centerStride 必须是正整数。');
    end
    if isempty(options.centerPhase)
        center = ceil((options.centerStride + 1) / 2);
        options.centerPhase = [center center];
    end
    if ~(isnumeric(options.centerPhase) && ...
            isequal(size(options.centerPhase), [1 2]) && ...
            all(isfinite(options.centerPhase)) && ...
            all(options.centerPhase == fix(options.centerPhase)) && ...
            all(options.centerPhase >= 1) && ...
            all(options.centerPhase <= options.centerStride))
        error('collocate_day:OperatorOptions', ...
            'centerPhase 必须是 1×2 整数，且位于 1 到 centerStride。');
    end
    if options.parameterWeighting ~= "none" && ~hasForwardModel
        error('collocate_day:OperatorOptions', ...
            '加权参数 TB 算子要求显式传入 MODEL。');
    end
end

function weights = equal_area_subgrid_weights(x, y, grid, divisions)
    fineCellSize = grid.cellSize / divisions;
    fineRows = grid.nRows * divisions;
    fineCols = grid.nCols * divisions;
    fineRow = floor((grid.yMax - y) ./ fineCellSize) + 1;
    fineCol = floor((x - grid.xMin) ./ fineCellSize) + 1;
    fineRow = min(max(fineRow, 1), fineRows);
    fineCol = min(max(fineCol, 1), fineCols);
    fineIndex = fineRow + (fineCol - 1) .* fineRows;
    [~, ~, group] = unique(fineIndex);
    samplesPerSubgrid = accumarray(group, 1);
    weights = 1 ./ samplesPerSubgrid(group);
end

function [weightedSums, denominator, support] = ...
        inverse_distance_center_accumulators( ...
        row, col, x, y, fields, grid, distancePower, windowSize)
%INVERSE_DISTANCE_CENTER_ACCUMULATORS 为每个窗口中心建立中心相关权重。
%   一个 OIB 点最多贡献给 windowSize^2 个候选 SMOS 中心；其权重由该点
%   到具体中心格几何中心的距离决定，不能复用 conv2 的固定点权重。

    shape = [grid.nRows, grid.nCols];
    nCells = prod(shape);
    nPoints = numel(row);
    candidatesPerPoint = windowSize ^ 2;
    targetIndex = zeros(candidatesPerPoint * nPoints, 1);
    sourceIndex = zeros(candidatesPerPoint * nPoints, 1);
    distanceM = zeros(candidatesPerPoint * nPoints, 1);
    cursor = 0;
    radius = (windowSize - 1) / 2;
    for rowOffset = -radius:radius
        centerRow = row - rowOffset;
        for colOffset = -radius:radius
            centerCol = col - colOffset;
            inside = centerRow >= 1 & centerRow <= grid.nRows & ...
                centerCol >= 1 & centerCol <= grid.nCols;
            pointIndex = find(inside);
            count = numel(pointIndex);
            if count == 0
                continue;
            end
            positions = cursor + (1:count);
            selectedRow = centerRow(inside);
            selectedCol = centerCol(inside);
            centerX = grid.xMin + ...
                (double(selectedCol) - 0.5) .* grid.cellSize;
            centerY = grid.yMax - ...
                (double(selectedRow) - 0.5) .* grid.cellSize;
            targetIndex(positions) = sub2ind( ...
                shape, selectedRow, selectedCol);
            sourceIndex(positions) = pointIndex;
            distanceM(positions) = hypot( ...
                x(pointIndex) - centerX, y(pointIndex) - centerY);
            cursor = cursor + count;
        end
    end
    targetIndex = targetIndex(1:cursor);
    sourceIndex = sourceIndex(1:cursor);
    distanceM = distanceM(1:cursor);

    zeroToleranceM = max(1e-6, 16 * eps(grid.cellSize));
    atCenter = distanceM <= zeroToleranceM;
    centerPointCount = accumarray( ...
        targetIndex, double(atCenter), [nCells 1], @sum, 0);
    hasCenterPoint = centerPointCount(targetIndex) > 0;
    weights = zeros(size(distanceM));
    regular = ~hasCenterPoint;
    weights(regular) = distanceM(regular) .^ (-distancePower);
    weights(hasCenterPoint & atCenter) = 1;

    denominatorVector = accumarray( ...
        targetIndex, weights, [nCells 1], @sum, 0);
    supportVector = accumarray( ...
        targetIndex, double(weights > 0), [nCells 1], @sum, 0);
    denominator = reshape(denominatorVector, shape);
    support = reshape(supportVector, shape);
    weightedSums = cell(1, numel(fields));
    for k = 1:numel(fields)
        value = accumarray(targetIndex, ...
            weights .* fields{k}(sourceIndex), [nCells 1], @sum, 0);
        weightedSums{k} = reshape(value, shape);
    end
end

function tb = evaluate_pair_forward( ...
        hi, hs, isMYI, isFYI, model, forwardFcn)
    tb = nan(size(hi));
    fyiModel = model;
    myiModel = model;
    if isfield(model, 'surfaceTempC') && ~isscalar(model.surfaceTempC)
        fyiModel.surfaceTempC = model.surfaceTempC(isFYI);
        myiModel.surfaceTempC = model.surfaceTempC(isMYI);
    end
    tb(isFYI) = forwardFcn(hi(isFYI), hs(isFYI), false, fyiModel);
    tb(isMYI) = forwardFcn(hi(isMYI), hs(isMYI), true, myiModel);
    if any(~isfinite(tb))
        error('collocate_day:ForwardTB', ...
            '配对尺度前向亮温包含非有限值。');
    end
end
