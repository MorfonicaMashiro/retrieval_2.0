function data = read_smap_l1c_product(filePath, productConfig, smapConfig)
%READ_SMAP_L1C_PRODUCT SMAP L1C 36/9 km 北极半轨的共享读取实现。

    filePath = char(string(filePath));
    if ~isfile(filePath)
        error('read_smap_l1c_product:MissingFile', ...
            'SMAP file does not exist: %s', filePath);
    end
    validate_reader_config(productConfig, smapConfig);
    [sourceFile, orbitNumber, orbitDirection, startTimeUtc] = ...
        parse_filename(filePath, productConfig);
    group = locate_north_group(filePath, productConfig.groupCandidates);

    latitude = read_numeric(filePath, group + "/cell_lat", 0, true);
    longitude = read_numeric(filePath, group + "/cell_lon", 0, true);
    n = numel(latitude);
    cellRow = read_optional_numeric( ...
        filePath, group + "/cell_row", n);
    if all(isnan(cellRow))
        cellRow = read_optional_numeric( ...
            filePath, group + "/cell_row_index", n);
    end
    cellColumn = read_optional_numeric( ...
        filePath, group + "/cell_column", n);
    if all(isnan(cellColumn))
        cellColumn = read_optional_numeric( ...
            filePath, group + "/cell_col", n);
    end
    require_length(longitude, n, 'cell_lon', filePath);
    gridSurfaceStatus = read_optional_numeric( ...
        filePath, group + "/cell_grid_surface_status", n);

    lookTables = cell(numel(smapConfig.looks), 1);
    for lookIndex = 1:numel(smapConfig.looks)
        look = string(smapConfig.looks(lookIndex));
        suffix = "_" + look;

        rawTbH = read_required_numeric( ...
            filePath, group + "/cell_tb_h" + suffix, n);
        rawTbV = read_required_numeric( ...
            filePath, group + "/cell_tb_v" + suffix, n);
        correctedTbH = read_optional_numeric( ...
            filePath, group + "/cell_tb_h_surface_corrected" + suffix, n);
        correctedTbV = read_optional_numeric( ...
            filePath, group + "/cell_tb_v_surface_corrected" + suffix, n);
        [tbH, tbV] = select_tb_fields(rawTbH, rawTbV, ...
            correctedTbH, correctedTbV, smapConfig.tbField, filePath);

        qualityFlagH = read_required_flag( ...
            filePath, group + "/cell_tb_qual_flag_h" + suffix, n);
        qualityFlagV = read_required_flag( ...
            filePath, group + "/cell_tb_qual_flag_v" + suffix, n);
        qualityH = required_bits_clear( ...
            qualityFlagH, smapConfig.requiredClearQualityBits);
        qualityV = required_bits_clear( ...
            qualityFlagV, smapConfig.requiredClearQualityBits);

        tbErrorH = read_optional_numeric( ...
            filePath, group + "/cell_tb_error_h" + suffix, n);
        tbErrorV = read_optional_numeric( ...
            filePath, group + "/cell_tb_error_v" + suffix, n);
        timeSecondsJ2000 = read_optional_numeric( ...
            filePath, group + "/cell_tb_time_seconds" + suffix, n);
        timeUtc = read_optional_time_utc( ...
            filePath, group + "/cell_tb_time_utc" + suffix, n);
        incidenceAngleDeg = read_optional_numeric( ...
            filePath, group + "/cell_boresight_incidence" + suffix, n);
        antennaScanAngleDeg = read_optional_numeric( ...
            filePath, group + "/cell_antenna_scan_angle" + suffix, n);
        centroidLatitude = read_optional_numeric( ...
            filePath, group + "/cell_lat_centroid" + suffix, n);
        centroidLongitude = read_optional_numeric( ...
            filePath, group + "/cell_lon_centroid" + suffix, n);
        nMeasurementsH = read_optional_numeric( ...
            filePath, group + "/cell_number_measurements_h" + suffix, n);
        nMeasurementsV = read_optional_numeric( ...
            filePath, group + "/cell_number_measurements_v" + suffix, n);
        waterFractionH = read_optional_numeric(filePath, ...
            group + "/cell_surface_water_fraction_mb_h" + suffix, n);
        waterFractionV = read_optional_numeric(filePath, ...
            group + "/cell_surface_water_fraction_mb_v" + suffix, n);

        validGeolocation = isfinite(latitude) & isfinite(longitude) & ...
            latitude >= 0 & latitude <= 90 & ...
            longitude >= -180 & longitude <= 180;
        rangeK = double(smapConfig.validTBRangeK);
        validH = validGeolocation & isfinite(tbH) & ...
            tbH >= rangeK(1) & tbH <= rangeK(2) & qualityH;
        validV = validGeolocation & isfinite(tbV) & ...
            tbV >= rangeK(1) & tbV <= rangeK(2) & qualityV;
        validBoth = validH & validV;
        if smapConfig.requireBothPolarizations
            recommendedMask = validBoth;
        else
            recommendedMask = validH | validV;
        end

        lookColumn = repmat(look, n, 1);
        sourceColumn = repmat(sourceFile, n, 1);
        productColumn = repmat(string(productConfig.product), n, 1);
        resolutionKm = repmat(double(productConfig.resolutionKm), n, 1);
        lookTables{lookIndex} = table(latitude, longitude, cellRow, ...
            cellColumn, gridSurfaceStatus, lookColumn, rawTbH, rawTbV, ...
            correctedTbH, correctedTbV, tbH, tbV, tbErrorH, tbErrorV, ...
            qualityFlagH, qualityFlagV, timeSecondsJ2000, timeUtc, ...
            incidenceAngleDeg, antennaScanAngleDeg, centroidLatitude, ...
            centroidLongitude, nMeasurementsH, nMeasurementsV, ...
            waterFractionH, waterFractionV, validH, validV, validBoth, ...
            recommendedMask, sourceColumn, productColumn, resolutionKm, ...
            'VariableNames', {'latitude','longitude','cellRow', ...
            'cellColumn','gridSurfaceStatus','look','rawTbH','rawTbV', ...
            'surfaceCorrectedTbH','surfaceCorrectedTbV','tbH','tbV', ...
            'tbErrorH','tbErrorV','qualityFlagH','qualityFlagV', ...
            'timeSecondsJ2000','timeUtc','incidenceAngleDeg', ...
            'antennaScanAngleDeg','centroidLatitude','centroidLongitude', ...
            'nMeasurementsH','nMeasurementsV','waterFractionH', ...
            'waterFractionV','validH','validV','validBoth', ...
            'recommendedMask','sourceFile','product','resolutionKm'});
    end

    observations = vertcat(lookTables{:});
    data = struct();
    data.source = string(filePath);
    data.sourceFile = sourceFile;
    data.product = string(productConfig.product);
    data.shortName = string(productConfig.shortName);
    data.version = string(productConfig.version);
    data.resolutionKm = double(productConfig.resolutionKm);
    data.grid = string(productConfig.grid);
    data.northEpsg = double(productConfig.northEpsg);
    data.northGroup = group;
    data.orbitNumber = orbitNumber;
    data.orbitDirection = orbitDirection;
    data.startTimeUtc = startTimeUtc;
    data.day = dateshift(startTimeUtc, 'start', 'day');
    data.tbField = string(smapConfig.tbField);
    data.observations = observations;
    data.recommended = observations(observations.recommendedMask, :);
end

function validate_reader_config(productConfig, smapConfig)
    productFields = {'product','shortName','version','resolutionKm', ...
        'grid','northEpsg','filePrefix','groupCandidates'};
    if ~all(isfield(productConfig, productFields))
        error('read_smap_l1c_product:ProductConfig', ...
            'SMAP product configuration is incomplete.');
    end
    if ~ismember(string(smapConfig.tbField), ...
            ["uncorrected", "surface_corrected"])
        error('read_smap_l1c_product:TBField', ...
            'smap.tbField must be uncorrected or surface_corrected.');
    end
    if ~isnumeric(smapConfig.validTBRangeK) || ...
            ~isequal(size(smapConfig.validTBRangeK), [1 2]) || ...
            smapConfig.validTBRangeK(2) <= smapConfig.validTBRangeK(1)
        error('read_smap_l1c_product:TBRange', ...
            'smap.validTBRangeK must be an increasing two-element vector.');
    end
end

function [sourceFile, orbitNumber, orbitDirection, startTimeUtc] = ...
        parse_filename(filePath, productConfig)
    [~, stem, extension] = fileparts(filePath);
    sourceFile = string(stem) + string(extension);
    expression = "^" + string(productConfig.filePrefix) + ...
        "_(\d{5})_([AD])_(\d{8}T\d{6})_R\d+_\d+\.h5$";
    tokens = regexp(sourceFile, expression, 'tokens', 'once');
    if isempty(tokens)
        error('read_smap_l1c_product:Filename', ...
            'Unexpected %s filename: %s', ...
            string(productConfig.product), sourceFile);
    end
    orbitNumber = str2double(tokens{1});
    orbitDirection = string(tokens{2});
    startTimeUtc = datetime(tokens{3}, ...
        'InputFormat', "yyyyMMdd'T'HHmmss", 'TimeZone', 'UTC');
end

function group = locate_north_group(filePath, candidates)
    info = h5info(filePath);
    groupNames = string({info.Groups.Name});
    candidates = string(candidates);
    selected = candidates(ismember(candidates, groupNames));
    if isempty(selected)
        error('read_smap_l1c_product:NorthGroup', ...
            'No north polar group in %s. Groups: %s', ...
            filePath, strjoin(groupNames, ', '));
    end
    group = selected(1);
end

function values = read_required_numeric(filePath, datasetPath, n)
    values = read_numeric(filePath, datasetPath, n, true);
    require_length(values, n, datasetPath, filePath);
end

function values = read_optional_numeric(filePath, datasetPath, n)
    values = read_numeric(filePath, datasetPath, n, false);
    values = normalize_optional_length(values, n, datasetPath, filePath);
end

function values = read_numeric(filePath, datasetPath, n, required)
    if ~dataset_exists(filePath, datasetPath)
        if required
            error('read_smap_l1c_product:MissingDataset', ...
                'Missing SMAP dataset %s in %s.', datasetPath, filePath);
        end
        values = nan(n, 1);
        return
    end
    values = double(h5read(filePath, char(datasetPath)));
    values = values(:);
    values(values <= -9990) = NaN;
end

function flags = read_required_flag(filePath, datasetPath, n)
    if ~dataset_exists(filePath, datasetPath)
        error('read_smap_l1c_product:MissingDataset', ...
            'Missing SMAP quality flag %s in %s.', datasetPath, filePath);
    end
    flags = uint16(h5read(filePath, char(datasetPath)));
    flags = flags(:);
    require_length(flags, n, datasetPath, filePath);
end

function timeUtc = read_optional_time_utc(filePath, datasetPath, n)
    timeUtc = NaT(n, 1, 'TimeZone', 'UTC');
    if ~dataset_exists(filePath, datasetPath)
        return
    end
    text = string(h5read(filePath, char(datasetPath)));
    text = text(:);
    require_length(text, n, datasetPath, filePath);
    valid = strlength(text) > 0 & ~contains(text, "-9999");
    if any(valid)
        try
            timeUtc(valid) = datetime(text(valid), ...
                'InputFormat', "yyyy-MM-dd'T'HH:mm:ss.SSS'Z'", ...
                'TimeZone', 'UTC');
        catch errorInfo
            error('read_smap_l1c_product:TimeFormat', ...
                'Invalid UTC text in %s: %s', datasetPath, errorInfo.message);
        end
    end
end

function clear = required_bits_clear(flags, zeroBasedBits)
    clear = true(size(flags));
    for bitIndex = reshape(double(zeroBasedBits), 1, [])
        clear = clear & bitget(flags, bitIndex + 1) == 0;
    end
end

function [tbH, tbV] = select_tb_fields(rawH, rawV, correctedH, correctedV, ...
        fieldName, filePath)
    fieldName = string(fieldName);
    if fieldName == "uncorrected"
        tbH = rawH;
        tbV = rawV;
        return
    end
    if all(isnan(correctedH)) || all(isnan(correctedV))
        error('read_smap_l1c_product:MissingSurfaceCorrectedTB', ...
            'Surface-corrected TB was requested but is absent in %s.', filePath);
    end
    tbH = correctedH;
    tbV = correctedV;
end

function values = normalize_optional_length(values, n, fieldName, filePath)
    if isempty(values) && n == 0
        values = nan(0, 1);
        return
    end
    require_length(values, n, fieldName, filePath);
end

function require_length(values, n, fieldName, filePath)
    if numel(values) ~= n
        error('read_smap_l1c_product:FieldSize', ...
            'SMAP field %s has %d values; expected %d in %s.', ...
            fieldName, numel(values), n, filePath);
    end
end

function present = dataset_exists(filePath, datasetPath)
    try
        h5info(filePath, char(datasetPath));
        present = true;
    catch errorInfo
        if contains(errorInfo.identifier, 'h5info:unableToFind') || ...
                contains(errorInfo.identifier, 'hdf5lib2') || ...
                contains(errorInfo.message, 'not found') || ...
                contains(errorInfo.message, 'does not exist')
            present = false;
        else
            rethrow(errorInfo)
        end
    end
end
