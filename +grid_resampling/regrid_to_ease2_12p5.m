function result = regrid_to_ease2_12p5(data, sourceKind, options)
%REGRID_TO_EASE2_12P5 统一进入 EASE2-N 12.5 km 空间索引体系。
%   RESULT = GRID_RESAMPLING.REGRID_TO_EASE2_12P5(DATA) 自动识别
%   READ_OIB、READ_SMOS、READ_SMAP_36KM 或 READ_SMAP_9KM 的输出。
%
%   SOURCEKIND 可显式设为 "oib"、"smos" 或 "smap"。目标网格固定为
%   EASE2-N 12.5 km：EPSG:6931、1440×1440、格心间距 12,500 m。
%
%   三类数据采用不同处理：
%     OIB  航迹点：保留每条原始记录，仅附加目标行列号与投影坐标；
%     SMAP 观测点：同样保留逐观测记录、前后视和 UTC，不做时间平均；
%     SMOS 规则栅格：按目标格心执行最近邻源像元映射，返回目标栅格。
%
%   OPTIONS.SMAPSELECTION="all"（默认）保留全部 SMAP 记录；设为
%   "recommended" 时使用读取器提供的推荐质量子集。

    if nargin < 2 || isempty(sourceKind)
        sourceKind = "auto";
    end
    if nargin < 3 || isempty(options)
        options = struct();
    end

    % 先确定传感器类型并补齐可选参数，避免在各分支内重复判断。
    sourceKind = normalize_source_kind(data, sourceKind);
    options = validate_options(options);
    grid = grid_resampling.ease_grid("EASE2");

    % 点数据只定位、不聚合；栅格数据调用专用最近邻映射。
    switch sourceKind
        case "oib"
            result = locate_oib_records(data, grid);

        case "smos"
            result = grid_resampling.regrid_smos_to_ease( ...
                data, "EASE2");
            result.sourceKind = "smos";
            result.targetGridName = "EASE2-N12.5km";
            result.targetSpacingKm = 12.5;
            result.physicalResolutionChanged = false;

        case "smap"
            result = locate_smap_records( ...
                data, grid, options.smapSelection);

        otherwise
            % normalize_source_kind 正常情况下已拦截该分支。
            error('regrid_to_ease2_12p5:SourceKind', ...
                'sourceKind must be oib, smos, or smap.');
    end
end

function sourceKind = normalize_source_kind(data, sourceKind)
%NORMALIZE_SOURCE_KIND 校验显式类型，或依据读取器输出结构自动识别。

    sourceKind = lower(string(sourceKind));
    if ~isscalar(sourceKind) || ...
            ~ismember(sourceKind, ["auto", "oib", "smos", "smap"])
        error('regrid_to_ease2_12p5:SourceKind', ...
            'sourceKind must be auto, oib, smos, or smap.');
    end
    if sourceKind ~= "auto"
        return
    end

    % 自动识别只依赖各读取器稳定的顶层接口，不依赖文件名：
    % OIB 有 table；SMOS 有二维 tb/latitude/longitude/nPair；
    % SMAP 有 observations 表。
    if isstruct(data) && isfield(data, 'table') && istable(data.table)
        sourceKind = "oib";
    elseif isstruct(data) && ...
            all(isfield(data, {'tb','latitude','longitude','nPair'}))
        sourceKind = "smos";
    elseif isstruct(data) && isfield(data, 'observations') && ...
            istable(data.observations)
        sourceKind = "smap";
    else
        error('regrid_to_ease2_12p5:AutoDetect', ...
            'Could not identify data as OIB, SMOS, or SMAP reader output.');
    end
end

function options = validate_options(options)
%VALIDATE_OPTIONS 只接受已声明选项，防止拼写错误被静默忽略。

    defaults = struct('smapSelection', "all");
    if ~isstruct(options) || ~isscalar(options)
        error('regrid_to_ease2_12p5:Options', ...
            'options must be a scalar struct.');
    end
    names = fieldnames(options);
    unknown = setdiff(names, fieldnames(defaults));
    if ~isempty(unknown)
        error('regrid_to_ease2_12p5:Options', ...
            'Unknown option: %s.', strjoin(unknown, ', '));
    end

    % 把用户给出的字段覆盖到默认结构。当前只有 smapSelection。
    for index = 1:numel(names)
        defaults.(names{index}) = options.(names{index});
    end
    options = defaults;
    options.smapSelection = lower(string(options.smapSelection));
    if ~isscalar(options.smapSelection) || ...
            ~ismember(options.smapSelection, ["all", "recommended"])
        error('regrid_to_ease2_12p5:SmapSelection', ...
            'options.smapSelection must be all or recommended.');
    end
end

function result = locate_oib_records(data, grid)
%LOCATE_OIB_RECORDS 为每条 OIB 航迹记录附加 EASE2 行列与投影坐标。

    if ~isstruct(data) || ~isfield(data, 'table') || ...
            ~istable(data.table)
        error('regrid_to_ease2_12p5:OIB', ...
            'OIB input must be a struct returned by read_oib.');
    end
    records = data.table;
    required = {'lat','lon'};
    if ~all(ismember(required, records.Properties.VariableNames))
        error('regrid_to_ease2_12p5:OIB', ...
            'OIB table must contain lat and lon.');
    end

    % 不平均 OIB 物理量；只定位每条记录，并另外统计每格点数。
    [records, count] = locate_table(records, 'lat', 'lon', grid);
    source = "";
    if isfield(data, 'source')
        source = string(data.source);
    end
    result = common_point_result(records, count, grid, "oib");
    result.source = source;
    result.sourceResolutionKm = NaN;
    result.sourceSampling = "airborne_track_points";
    result.method = "point_to_cell_without_value_aggregation";
end

function result = locate_smap_records(data, grid, selection)
%LOCATE_SMAP_RECORDS 定位 SMAP 逐观测记录，不提前合成前视、后视或时间。

    if ~isstruct(data) || ~isfield(data, 'observations') || ...
            ~istable(data.observations)
        error('regrid_to_ease2_12p5:SMAP', ...
            'SMAP input must be returned by read_smap_36km or read_smap_9km.');
    end

    % recommended 由 SMAP 读取器根据产品质量标志生成；all 保留全部记录。
    if selection == "recommended"
        if ~isfield(data, 'recommended') || ~istable(data.recommended)
            error('regrid_to_ease2_12p5:SMAP', ...
                'SMAP input does not contain a recommended table.');
        end
        records = data.recommended;
    else
        records = data.observations;
    end
    required = {'latitude','longitude'};
    if ~all(ismember(required, records.Properties.VariableNames))
        error('regrid_to_ease2_12p5:SMAP', ...
            'SMAP observations must contain latitude and longitude.');
    end

    [records, count] = locate_table( ...
        records, 'latitude', 'longitude', grid);
    result = common_point_result(records, count, grid, "smap");
    result.selection = selection;
    result.source = string(data.source);
    result.product = string(data.product);
    result.sourceResolutionKm = double(data.resolutionKm);
    result.sourceSampling = "satellite_footprints";
    result.method = "point_to_cell_without_time_aggregation";
    result.resolutionNote = ...
        "12.5 km is target spacing; native physical resolution is retained in metadata.";
end

function [records, count] = locate_table( ...
        records, latitudeName, longitudeName, grid)
%LOCATE_TABLE 定位表中经纬度，并生成每个 EASE2 格的记录数。

    latitude = double(records.(latitudeName));
    longitude = double(records.(longitudeName));
    [row, col, x, y] = grid_resampling.PSP2EASE( ...
        latitude, longitude, "EASE2");

    % 统一成列向量，以满足 addvars 对表高的要求；原记录顺序不改变。
    row = double(row(:));
    col = double(col(:));
    x = double(x(:));
    y = double(y(:));
    records = addvars(records, row, col, x, y, ...
        'NewVariableNames', ...
        {'ease2Row','ease2Col','ease2X','ease2Y'});

    % 仅对落入目标网格的记录建立线性索引，并按像元累加记录数。
    % count 是覆盖诊断，不是观测值平均结果。
    inside = isfinite(row) & isfinite(col);
    countVector = zeros(grid.nRows * grid.nCols, 1, 'uint32');
    if any(inside)
        cellIndex = sub2ind( ...
            [grid.nRows, grid.nCols], row(inside), col(inside));
        countVector = uint32(accumarray( ...
            cellIndex, 1, [numel(countVector), 1], @sum, 0));
    end
    count = reshape(countVector, grid.nRows, grid.nCols);
end

function result = common_point_result(records, count, grid, sourceKind)
%COMMON_POINT_RESULT 组装 OIB/SMAP 点定位结果的统一输出结构。

    inside = isfinite(records.ease2Row) & ...
        isfinite(records.ease2Col);
    result = struct();
    result.sourceKind = sourceKind;
    result.targetGridName = "EASE2-N12.5km";
    result.targetSpacingKm = 12.5;
    result.grid = grid;

    % records 保留网格外记录用于审计；validRecords 仅保留成功定位记录。
    result.records = records;
    result.validRecords = records(inside, :);
    result.count = count;

    % 点落格不改变传感器足迹或航迹采样分辨率。
    result.physicalResolutionChanged = false;
end
