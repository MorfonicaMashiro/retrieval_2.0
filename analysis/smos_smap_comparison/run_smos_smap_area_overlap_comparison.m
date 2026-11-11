function [pairs, metrics, summary, overlapWeights] = ...
        run_smos_smap_area_overlap_comparison()
%RUN_SMOS_SMAP_AREA_OVERLAP_COMPARISON 两日 36 km 共同支撑域比较。
%   每个 SMAP SPL1CTB 36 km EASE2 网格是目标空间支撑域；SMOS 先重建
%   到原版 EASE 12.53376 km，再按真实多边形交叠面积汇总到该目标格。
%   SMAP 同时输出最接近逐格 OIB 中位时刻的单次观测和同日全部观测均值。

    branchRoot = fileparts(mfilename('fullpath'));
    matlabRoot = fileparts(fileparts(branchRoot));
    workspaceRoot = fileparts(matlabRoot);
    dataRoot = fullfile(workspaceRoot, 'data');
    resultRoot = fullfile(branchRoot, 'results');
    addpath(genpath(matlabRoot));
    if ~isfolder(resultRoot), mkdir(resultRoot); end

    config = retrieval_config();
    targetDays = [datetime(2015,4,1); datetime(2015,4,3)];
    outputStem = "smos_smap_20150401_20150403_smap36_area_overlap";
    minimumSmosCoverage = ...
        config.data.smap.comparison.minimumSmosAreaCoverageFraction;
    operatorOptions = struct('centerStride', 3, 'centerPhase', [2 2]);
    overlapOptions = struct('edgeSegments', ...
        config.data.smap.comparison.overlapBoundarySegments, ...
        'areaClosureTolerance', 2e-5);

    pairParts = cell(numel(targetDays), 1);
    weightParts = cell(numel(targetDays), 1);
    summaryParts = cell(numel(targetDays), 1);
    inventoryParts = cell(numel(targetDays), 1);
    sampleParts = cell(numel(targetDays), 1);

    for dayIndex = 1:numel(targetDays)
        currentDay = targetDays(dayIndex);
        token = char(string(currentDay, 'yyyyMMdd'));
        oibFile = fullfile(dataRoot, 'OIB', 'NSIDC-0708_QuickLook', ...
            '2015_GR_NASA', ['OIB_' token '_IDCSI2.txt']);
        smosFile = fullfile(dataRoot, 'smos', '2015', ...
            ['SMOS_TB_north_' token '.nc']);
        smapDirectory = fullfile(dataRoot, 'smap', 'SPL1CTB_006', ...
            '2015', char(string(currentDay,'MM')), ...
            char(string(currentDay,'dd')));
        require_file(oibFile, 'OIB');
        require_file(smosFile, 'SMOS');
        smapListing = dir(fullfile(smapDirectory, ...
            ['SMAP_L1C_TB_*_' token 'T*_R19240_001.h5']));
        if isempty(smapListing)
            error('run_smos_smap_area_overlap_comparison:MissingSMAP', ...
                '%s has no SMAP SPL1CTB files in %s.', token, smapDirectory);
        end
        smapFiles = sort(string(fullfile( ...
            {smapListing.folder}, {smapListing.name})).');

        oib = read_oib(oibFile);
        smosEase = grid_resampling.regrid_smos_to_ease( ...
            read_smos(smosFile), "EASE");
        [stridePairs, collocateStats] = collocate_day( ...
            oib, smosEase, currentDay, [], [], operatorOptions);
        targets = build_target_cells(oib, currentDay);
        [geometryWeights, geometryDiagnostics] = ...
            compute_smap36_ease_overlap_weights( ...
            targets.smapRow, targets.smapCol, overlapOptions);
        [targets, dayWeights] = aggregate_smos_to_targets( ...
            targets, geometryWeights, geometryDiagnostics, ...
            smosEase, currentDay, minimumSmosCoverage);
        [smapObservations, inventory] = read_smap_day( ...
            smapFiles, config.data.smap);
        targets = attach_smap_time_products(targets, smapObservations, ...
            config.data.smap.primaryTimeWindowHours, ...
            config.data.smap.sensitivityTimeWindowHours);

        targets.day = repmat(currentDay, height(targets), 1);
        targets.TBSMAPNearestMeanHV = ...
            (targets.TBSMAPNearestH + targets.TBSMAPNearestV) ./ 2;
        targets.TBSMAPDailyMeanHV = ...
            (targets.TBSMAPDailyH + targets.TBSMAPDailyV) ./ 2;
        targets.nearestMinusSMOS = ...
            targets.TBSMAPNearestMeanHV - targets.TBSMOSAreaWeighted;
        targets.dailyMinusSMOS = ...
            targets.TBSMAPDailyMeanHV - targets.TBSMOSAreaWeighted;
        pairParts{dayIndex} = targets;
        weightParts{dayIndex} = dayWeights;
        inventoryParts{dayIndex} = inventory;
        sampleParts{dayIndex} = observation_sample( ...
            smapObservations, targets, currentDay);

        adequate = targets.smosCoverageFraction >= minimumSmosCoverage;
        nearestAvailable = adequate & ...
            isfinite(targets.TBSMAPNearestMeanHV);
        within3 = nearestAvailable & targets.smapWithin3Hours;
        within6 = nearestAvailable & targets.smapWithin6Hours;
        summaryParts{dayIndex} = table(currentDay, numel(smapFiles), ...
            collocateStats.rawPoints, collocateStats.qualityPoints, ...
            height(stridePairs), height(targets), nnz(adequate), ...
            nnz(nearestAvailable), nnz(within3), nnz(within6), ...
            targets.flightStartUtc(1), targets.flightEndUtc(1), ...
            min(targets.oibStartUtc), max(targets.oibEndUtc), ...
            median(targets.smapNearestAbsDeltaHours(nearestAvailable), ...
            'omitnan'), median(targets.smosCoverageFraction, 'omitnan'), ...
            mean(targets.nearestMinusSMOS(within3), 'omitnan'), ...
            mean(targets.dailyMinusSMOS(adequate), 'omitnan'), ...
            'VariableNames', {'day','nSMAPFiles','rawOIBPoints', ...
            'qualityOIBPoints','legacyStride3Centers','targetSMAPCells', ...
            'adequateSMOSCoverageCells','nearestSMAPCells', ...
            'nearestWithin3Hours','nearestWithin6Hours', ...
            'flightStartUtc','flightEndUtc','targetOIBStartUtc', ...
            'targetOIBEndUtc','medianNearestAbsDeltaHours', ...
            'medianSMOSCoverageFraction','biasNearest3hSMAPMinusSMOS', ...
            'biasDailySMAPMinusSMOS'});
        fprintf(['%s: target SMAP cells=%d, SMOS coverage>=%.0f%%: %d, ' ...
            'nearest SMAP within 3 h: %d.\n'], token, height(targets), ...
            minimumSmosCoverage * 100, nnz(adequate), nnz(within3));
    end

    pairs = sortrows(vertcat(pairParts{:}), {'day','smapRow','smapCol'});
    overlapWeights = sortrows(vertcat(weightParts{:}), ...
        {'day','smapRow','smapCol','easeRow','easeCol'});
    summary = vertcat(summaryParts{:});
    fileInventory = vertcat(inventoryParts{:});
    smapObservationSample = vertcat(sampleParts{:});
    metrics = build_metrics(pairs, minimumSmosCoverage);
    method = struct( ...
        'targetSupport', "SMAP SPL1CTB.006 EASE2-N 36 km cell", ...
        'smosSourceGrid', "original EASE-N 12.53376 km", ...
        'overlapMethod', "densified polygon intersection in EPSG:6931", ...
        'edgeSegments', overlapOptions.edgeSegments, ...
        'minimumSmosCoverage', minimumSmosCoverage, ...
        'primarySmapTime', "single observation nearest OIB cell median UTC", ...
        'smapSensitivityTime', "mean of all valid same-day observations", ...
        'smosTime', "daily L3B composite; no pixel-level UTC");

    matPath = fullfile(resultRoot, outputStem + ".mat");
    save(matPath, 'pairs', 'metrics', 'summary', 'overlapWeights', ...
        'fileInventory', 'smapObservationSample', 'method', ...
        'operatorOptions', '-v7.3');
    writetable(pairs, fullfile(resultRoot, outputStem + "_pairs.csv"));
    writetable(metrics, fullfile(resultRoot, outputStem + "_metrics.csv"));
    writetable(summary, fullfile(resultRoot, outputStem + "_summary.csv"));
    writetable(overlapWeights, ...
        fullfile(resultRoot, outputStem + "_overlap_weights.csv"));
    writetable(fileInventory, ...
        fullfile(resultRoot, outputStem + "_smap_file_inventory.csv"));
    writetable(smapObservationSample, ...
        fullfile(resultRoot, outputStem + "_smap_observation_sample.csv"));
    fprintf('Area-overlap SMOS-SMAP comparison complete: %s\n', matPath);
end

function targets = build_target_cells(oib, currentDay)
    smapProjection = projcrs(6931);
    tbl = oib.table;
    elapsedAll = double(tbl.elapsed);
    dayUtc = datetime(string(currentDay, 'yyyy-MM-dd'), ...
        'InputFormat', 'yyyy-MM-dd', 'TimeZone', 'UTC');
    flightTime = dayUtc + seconds(elapsedAll(isfinite(elapsedAll)));
    flightStartUtc = min(flightTime);
    flightEndUtc = max(flightTime);

    fbSnow = double(tbl.mean_fb);
    snowDepth = double(tbl.snow_depth);
    thickness = double(tbl.thickness);
    iceType = double(tbl.my_ice_flag);
    valid = isfinite(tbl.lat) & isfinite(tbl.lon) & ...
        isfinite(elapsedAll) & isfinite(fbSnow) & ...
        isfinite(snowDepth) & isfinite(thickness) & ...
        isfinite(fbSnow - snowDepth) & tbl.pcnt_ow == 0 & ...
        tbl.pcnt_thin_ice == 0 & tbl.pcnt_grey_ice == 0 & ...
        ismember(iceType, [0 1]);
    longitude = mod(double(tbl.lon(valid)) + 180, 360) - 180;
    [pointX, pointY] = projfwd(smapProjection, ...
        double(tbl.lat(valid)), longitude);
    pointRow = floor((9000000 - pointY) ./ 36000) + 1;
    pointCol = floor((pointX + 9000000) ./ 36000) + 1;
    located = pointRow >= 1 & pointRow <= 500 & ...
        pointCol >= 1 & pointCol <= 500 & ...
        isfinite(pointRow) & isfinite(pointCol);
    pointRow = pointRow(located);
    pointCol = pointCol(located);
    pointTime = dayUtc + seconds(elapsedAll(valid));
    pointTime = pointTime(located);
    pointIceType = iceType(valid);
    pointIceType = pointIceType(located);
    uniqueCells = unique([pointRow pointCol], 'rows', 'sorted');
    n = size(uniqueCells, 1);
    targets = table(uniqueCells(:,1), uniqueCells(:,2), ...
        'VariableNames', {'smapRow','smapCol'});
    oibStartUtc = NaT(n, 1, 'TimeZone', 'UTC');
    oibMedianUtc = NaT(n, 1, 'TimeZone', 'UTC');
    oibEndUtc = NaT(n, 1, 'TimeZone', 'UTC');
    nOIBTimePoints = zeros(n, 1);
    nFYI = zeros(n,1);
    nMYI = zeros(n,1);
    for index = 1:n
        selected = pointRow == targets.smapRow(index) & ...
            pointCol == targets.smapCol(index);
        nOIBTimePoints(index) = nnz(selected);
        nFYI(index) = nnz(pointIceType(selected) == 0);
        nMYI(index) = nnz(pointIceType(selected) == 1);
        if any(selected)
            oibStartUtc(index) = min(pointTime(selected));
            oibMedianUtc(index) = median(pointTime(selected));
            oibEndUtc(index) = max(pointTime(selected));
        end
    end
    targets = addvars(targets, nOIBTimePoints, nFYI, nMYI, ...
        nFYI > 0 & nMYI > 0, oibStartUtc, oibMedianUtc, oibEndUtc, ...
        repmat(flightStartUtc,n,1), repmat(flightEndUtc,n,1), ...
        'NewVariableNames', {'nOIBTimePoints','nFYI','nMYI', ...
        'mixedIceType','oibStartUtc','oibMedianUtc','oibEndUtc', ...
        'flightStartUtc','flightEndUtc'});
end

function [targets, dayWeights] = aggregate_smos_to_targets( ...
        targets, weights, diagnostics, smosEase, currentDay, minimumCoverage)
    n = height(targets);
    tb = nan(n,1);
    coverage = zeros(n,1);
    nOverlap = zeros(n,1);
    nValid = zeros(n,1);
    geometryCoverage = zeros(n,1);
    dayParts = cell(n,1);
    targetArea = 36000 ^ 2;
    for index = 1:n
        selected = weights.smapRow == targets.smapRow(index) & ...
            weights.smapCol == targets.smapCol(index);
        current = weights(selected,:);
        linear = sub2ind(size(smosEase.tb), ...
            current.easeRow, current.easeCol);
        current.smosTB = smosEase.tb(linear);
        current.smosValid = isfinite(current.smosTB);
        validArea = sum(current.overlapAreaM2(current.smosValid));
        coverage(index) = validArea / targetArea;
        nOverlap(index) = height(current);
        nValid(index) = nnz(current.smosValid);
        if validArea > 0
            current.normalizedValidWeight = zeros(height(current),1);
            current.normalizedValidWeight(current.smosValid) = ...
                current.overlapAreaM2(current.smosValid) ./ validArea;
            tb(index) = sum(current.smosTB(current.smosValid) .* ...
                current.normalizedValidWeight(current.smosValid));
        else
            current.normalizedValidWeight = zeros(height(current),1);
        end
        current.day = repmat(currentDay, height(current), 1);
        dayParts{index} = current;
        diagnosticRow = diagnostics.smapRow == targets.smapRow(index) & ...
            diagnostics.smapCol == targets.smapCol(index);
        geometryCoverage(index) = ...
            diagnostics.geometryCoverageFraction(diagnosticRow);
    end
    targets = addvars(targets, tb, coverage, geometryCoverage, ...
        nOverlap, nValid, coverage >= minimumCoverage, ...
        'NewVariableNames', {'TBSMOSAreaWeighted', ...
        'smosCoverageFraction','geometryCoverageFraction', ...
        'nOverlapEASECells','nValidSMOSEASECells','adequateSMOSCoverage'});
    dayWeights = vertcat(dayParts{:});
end

function [observations, inventory] = read_smap_day(filePaths, smapConfig)
    observationParts = cell(numel(filePaths),1);
    inventoryParts = cell(numel(filePaths),1);
    projection = projcrs(6931);
    for index = 1:numel(filePaths)
        data = read_smap_36km(filePaths(index), smapConfig);
        current = data.recommended;
        current.orbitNumber = repmat(data.orbitNumber, height(current), 1);
        current.orbitDirection = repmat( ...
            data.orbitDirection, height(current), 1);
        current.fileStartUtc = repmat( ...
            data.startTimeUtc, height(current), 1);
        current.smapRow = double(current.cellRow) + 1;
        current.smapCol = double(current.cellColumn) + 1;
        [x, y] = projfwd(projection, current.latitude, current.longitude);
        calculatedRow = floor((9000000 - y) ./ 36000) + 1;
        calculatedCol = floor((x + 9000000) ./ 36000) + 1;
        if any(calculatedRow ~= current.smapRow | ...
                calculatedCol ~= current.smapCol)
            error('run_smos_smap_area_overlap_comparison:SMAPIndex', ...
                'Stored SMAP cell indices disagree with EPSG:6931 geometry.');
        end
        observationParts{index} = current;
        finiteTime = ~isnat(current.timeUtc);
        if any(finiteTime)
            firstUtc = min(current.timeUtc(finiteTime));
            lastUtc = max(current.timeUtc(finiteTime));
        else
            firstUtc = NaT(1,1,'TimeZone','UTC');
            lastUtc = NaT(1,1,'TimeZone','UTC');
        end
        inventoryParts{index} = table(data.sourceFile, data.orbitNumber, ...
            data.orbitDirection, data.startTimeUtc, ...
            height(data.observations), height(current), firstUtc, lastUtc, ...
            'VariableNames', {'sourceFile','orbitNumber','orbitDirection', ...
            'fileStartUtc','nObservations','nRecommended', ...
            'firstCellUtc','lastCellUtc'});
    end
    observations = vertcat(observationParts{:});
    inventory = vertcat(inventoryParts{:});
end

function targets = attach_smap_time_products( ...
        targets, observations, primaryHours, sensitivityHours)
    n = height(targets);
    nearestH = nan(n,1); nearestV = nan(n,1);
    dailyH = nan(n,1); dailyV = nan(n,1);
    nDaily = zeros(n,1);
    nearestUtc = NaT(n,1,'TimeZone','UTC');
    firstUtc = NaT(n,1,'TimeZone','UTC');
    lastUtc = NaT(n,1,'TimeZone','UTC');
    signedDelta = nan(n,1); absoluteDelta = nan(n,1);
    nearestLook = strings(n,1); nearestFile = strings(n,1);
    nearestOrbit = nan(n,1); nearestDirection = strings(n,1);
    nearestIncidence = nan(n,1);
    for index = 1:n
        selected = observations.smapRow == targets.smapRow(index) & ...
            observations.smapCol == targets.smapCol(index) & ...
            ~isnat(observations.timeUtc);
        current = observations(selected,:);
        nDaily(index) = height(current);
        if isempty(current), continue; end
        dailyH(index) = mean(current.tbH, 'omitnan');
        dailyV(index) = mean(current.tbV, 'omitnan');
        firstUtc(index) = min(current.timeUtc);
        lastUtc(index) = max(current.timeUtc);
        if isnat(targets.oibMedianUtc(index)), continue; end
        deltaHours = seconds(current.timeUtc - ...
            targets.oibMedianUtc(index)) ./ 3600;
        [absoluteDelta(index), nearestIndex] = min(abs(deltaHours));
        signedDelta(index) = deltaHours(nearestIndex);
        nearestH(index) = current.tbH(nearestIndex);
        nearestV(index) = current.tbV(nearestIndex);
        nearestUtc(index) = current.timeUtc(nearestIndex);
        nearestLook(index) = current.look(nearestIndex);
        nearestFile(index) = current.sourceFile(nearestIndex);
        nearestOrbit(index) = current.orbitNumber(nearestIndex);
        nearestDirection(index) = current.orbitDirection(nearestIndex);
        nearestIncidence(index) = current.incidenceAngleDeg(nearestIndex);
    end
    targets = addvars(targets, nDaily, firstUtc, lastUtc, dailyH, dailyV, ...
        nearestUtc, signedDelta, absoluteDelta, nearestH, nearestV, ...
        nearestLook, nearestFile, nearestOrbit, nearestDirection, ...
        nearestIncidence, absoluteDelta <= primaryHours, ...
        absoluteDelta <= sensitivityHours, ...
        'NewVariableNames', {'nSMAPObservationsDaily','smapFirstUtc', ...
        'smapLastUtc','TBSMAPDailyH','TBSMAPDailyV','smapNearestUtc', ...
        'smapNearestSignedDeltaHours','smapNearestAbsDeltaHours', ...
        'TBSMAPNearestH','TBSMAPNearestV','smapNearestLook', ...
        'smapNearestSourceFile','smapNearestOrbit','smapNearestDirection', ...
        'smapNearestIncidenceDeg','smapWithin3Hours','smapWithin6Hours'});
end

function metrics = build_metrics(pairs, minimumCoverage)
    periods = [NaT; unique(pairs.day)];
    strategies = ["nearest_all"; "nearest_3h"; ...
        "nearest_6h"; "same_day_mean"];
    channels = ["H"; "V"; "(H+V)/2"];
    metrics = table();
    for periodIndex = 1:numel(periods)
        if ismissing(periods(periodIndex))
            periodMask = true(height(pairs),1);
            periodLabel = "two_days";
        else
            periodMask = pairs.day == periods(periodIndex);
            periodLabel = string(periods(periodIndex), 'yyyy-MM-dd');
        end
        for strategyIndex = 1:numel(strategies)
            strategy = strategies(strategyIndex);
            switch strategy
                case "nearest_all"
                    strategyMask = true(height(pairs),1);
                    values = {pairs.TBSMAPNearestH, ...
                        pairs.TBSMAPNearestV, pairs.TBSMAPNearestMeanHV};
                case "nearest_3h"
                    strategyMask = pairs.smapWithin3Hours;
                    values = {pairs.TBSMAPNearestH, ...
                        pairs.TBSMAPNearestV, pairs.TBSMAPNearestMeanHV};
                case "nearest_6h"
                    strategyMask = pairs.smapWithin6Hours;
                    values = {pairs.TBSMAPNearestH, ...
                        pairs.TBSMAPNearestV, pairs.TBSMAPNearestMeanHV};
                otherwise
                    strategyMask = true(height(pairs),1);
                    values = {pairs.TBSMAPDailyH, pairs.TBSMAPDailyV, ...
                        pairs.TBSMAPDailyMeanHV};
            end
            for channelIndex = 1:numel(channels)
                selected = periodMask & strategyMask & ...
                    pairs.smosCoverageFraction >= minimumCoverage;
                x = pairs.TBSMOSAreaWeighted(selected);
                y = values{channelIndex}(selected);
                fit = fit_retrieval(x, y);
                metrics = [metrics; table(periodLabel, strategy, ...
                    channels(channelIndex), fit.n, fit.correlation, ...
                    fit.rSquared, fit.bias, fit.rmse, fit.slope, ...
                    fit.intercept, 'VariableNames', {'period', ...
                    'timeStrategy','SMAPChannel','n','correlation', ...
                    'rSquared','biasSMAPMinusSMOS','rmse','slope', ...
                    'intercept'})]; %#ok<AGROW>
            end
        end
    end
end

function sample = observation_sample(observations, targets, currentDay)
    targetKey = (targets.smapRow - 1) * 500 + targets.smapCol;
    observationKey = (observations.smapRow - 1) * 500 + ...
        observations.smapCol;
    selected = ismember(observationKey, targetKey);
    names = {'smapRow','smapCol','look','tbH','tbV','timeUtc', ...
        'incidenceAngleDeg','qualityFlagH','qualityFlagV','orbitNumber', ...
        'orbitDirection','sourceFile'};
    sample = observations(selected, names);
    sample = sortrows(sample, {'timeUtc','smapRow','smapCol'});
    sample = sample(1:min(20,height(sample)), :);
    sample.day = repmat(currentDay, height(sample), 1);
    sample = movevars(sample, 'day', 'Before', 1);
end

function require_file(filePath, label)
    if ~isfile(filePath)
        error('run_smos_smap_area_overlap_comparison:MissingInput', ...
            '%s file not found: %s', label, filePath);
    end
end
