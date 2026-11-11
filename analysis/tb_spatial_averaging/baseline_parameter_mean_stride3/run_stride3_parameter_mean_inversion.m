function [results, metricSummary, pairs, summary] = ...
        run_stride3_parameter_mean_inversion(targetDays, outputStem)
%RUN_STRIDE3_PARAMETER_MEAN_INVERSION 运行步长 3 的原始参数均值反演。
%   每个保留的 SMOS 中心格使用互不重叠的 3x3 EASE 窗口；窗口内先对
%   OIB hi/hs 作原始点算术平均，再调用 Burke 前向模型计算一次 TB。
%   TARGETDAYS 为空时处理本地 2012、2015 年全部 OIB 飞行日。

    if nargin < 1
        targetDays = datetime.empty(0, 1);
    end
    targetDays = dateshift(datetime(targetDays(:)), 'start', 'day');

    experimentRoot = fileparts(mfilename('fullpath'));
    matlabRoot = fileparts(fileparts(fileparts(experimentRoot)));
    dataRoot = fullfile(fileparts(matlabRoot), 'data');
    resultRoot = fullfile(experimentRoot, 'results');
    figureRoot = fullfile(experimentRoot, 'figures');
    addpath(matlabRoot, fullfile(matlabRoot, 'readdata'), ...
        fullfile(matlabRoot, 'plot'), experimentRoot);
    if ~isfolder(resultRoot), mkdir(resultRoot); end
    if ~isfolder(figureRoot), mkdir(figureRoot); end

    retrievalConfig = retrieval_config();
    retrievalConfig.collocation.computeSampleMeanTB = false;
    retrievalConfig.collocation.parameterWeighting = "none";
    model = retrievalConfig.model;
    operatorOptions = retrievalConfig.collocation;
    gridVersion = retrievalConfig.data.gridVersion;
    if nargin < 2 || strlength(string(outputStem)) == 0
        outputStem = "parameter_mean_" + ...
            retrievalConfig.output.variantTag + "_2012_2015";
    end
    outputStem = string(outputStem);
    flightFiles = find_flight_files(dataRoot, targetDays);
    if isempty(flightFiles)
        error('run_stride3_parameter_mean_inversion:NoFlights', ...
            '指定日期没有找到 2012/2015 OIB 飞行文件。');
    end

    pairs = table();
    summary = table();
    outputFile = fullfile(resultRoot, outputStem + ".mat");
    timer = tic;
    for flightIndex = 1:height(flightFiles)
        day = flightFiles.day(flightIndex);
        token = char(string(day, 'yyyyMMdd'));
        oibFile = char(flightFiles.path(flightIndex));
        smosFile = fullfile(dataRoot, 'smos', char(string(year(day))), ...
            ['SMOS_TB_north_' token '.nc']);
        if ~isfile(smosFile)
            error('run_stride3_parameter_mean_inversion:MissingSMOS', ...
                '%s 缺少同日 SMOS 文件：%s', token, smosFile);
        end

        oib = read_oib(oibFile);
        smos = read_smos(smosFile);
        smosEase = grid_resampling.regrid_smos_to_ease( ...
            smos, gridVersion);
        [dailyPairs, stats] = collocate_day( ...
            oib, smosEase, day, retrievalConfig);
        pairs = [pairs; dailyPairs]; %#ok<AGROW>
        summary = [summary; table(day, stats.rawPoints, ...
            stats.qualityPoints, stats.tbOverlapCells, stats.matchedCells, ...
            stats.mixedTypeCells, stats.centerStride, ...
            stats.centerPhase(1), stats.centerPhase(2), ...
            'VariableNames', {'day','rawPoints','qualityPoints', ...
            'tbOverlapCells','matchedCells','mixedTypeCells', ...
            'centerStride','centerPhaseRow','centerPhaseCol'})]; %#ok<AGROW>

        elapsedSeconds = toc(timer);
        etaSeconds = elapsedSeconds / flightIndex * ...
            (height(flightFiles) - flightIndex);
        fprintf(['[%d/%d, %5.1f%%] %s：步长 3 配对 %d 条；' ...
            '累计 %.1f min，预计剩余 %.1f min。\n'], ...
            flightIndex, height(flightFiles), ...
            100 * flightIndex / height(flightFiles), token, ...
            stats.matchedCells, elapsedSeconds / 60, etaSeconds / 60);
    end

    pairs = sortrows(pairs, {'day','row','col'});
    summary = sortrows(summary, 'day');
    referenceMatch = verify_against_stride_filter( ...
        pairs, matlabRoot, targetDays, operatorOptions);
    save(outputFile, 'pairs', 'summary', 'model', 'operatorOptions', ...
        'gridVersion', 'referenceMatch', 'retrievalConfig', '-v7.3');
    writetable(pairs, fullfile(resultRoot, outputStem + "_pairs.csv"));
    writetable(summary, fullfile(resultRoot, outputStem + "_summary.csv"));

    inputNames = {'day','flight','row','col','TB','FBice', ...
        'FBsnow','isMYI','hi','hs','nOIB'};
    if retrievalConfig.temperature.source == "oib_kt19"
        inputNames{end + 1} = 'surfaceTempC';
    end
    inputs = pairs(:, inputNames);
    tb_forward(1, 0.1, false, model);
    n = height(inputs);
    hiRadar = nan(n, 1); hsRadar = nan(n, 1);
    hiLaser = nan(n, 1); hsLaser = nan(n, 1);
    hiLaserAlt = nan(n, 1); hsLaserAlt = nan(n, 1);
    nRadarRoots = zeros(n, 1); nLaserRoots = zeros(n, 1);
    tbFromOIB = nan(n, 1);
    solverTimer = tic;
    progressEvery = max(1, ceil(n / 100));
    fprintf(['开始逐条反演：%d 条；Burke 递推，Vb>70‰ 标记但不截断，' ...
        '先平均 hi/hs 后计算 TB。\n'], n);
    for pairIndex = 1:n
        pairConfig = retrievalConfig;
        pairModel = model;
        if retrievalConfig.temperature.source == "oib_kt19"
            pairModel.surfaceTempC = inputs.surfaceTempC(pairIndex);
            pairConfig.model = pairModel;
        end
        [radar, ~] = invert_pair( ...
            inputs.TB(pairIndex), inputs.FBice(pairIndex), "ice", ...
            inputs.isMYI(pairIndex), pairConfig);
        [laser, ~] = invert_pair( ...
            inputs.TB(pairIndex), inputs.FBsnow(pairIndex), "snow", ...
            inputs.isMYI(pairIndex), pairConfig);
        nRadarRoots(pairIndex) = numel(radar);
        nLaserRoots(pairIndex) = numel(laser);
        if ~isempty(radar)
            hiRadar(pairIndex) = radar(1).hi;
            hsRadar(pairIndex) = radar(1).hs;
        end
        if ~isempty(laser)
            hiLaser(pairIndex) = laser(1).hi;
            hsLaser(pairIndex) = laser(1).hs;
        end
        if numel(laser) > 1
            hiLaserAlt(pairIndex) = laser(2).hi;
            hsLaserAlt(pairIndex) = laser(2).hs;
        end
        tbFromOIB(pairIndex) = tb_forward(inputs.hi(pairIndex), ...
            inputs.hs(pairIndex), inputs.isMYI(pairIndex), pairModel);

        if pairIndex == 1 || mod(pairIndex, progressEvery) == 0 || ...
                pairIndex == n
            elapsedSeconds = toc(solverTimer);
            etaSeconds = elapsedSeconds / pairIndex * (n - pairIndex);
            fprintf(['[%d/%d, %5.1f%%] %s；已用 %.1f min，' ...
                '预计剩余 %.1f min；FBice 有解 %d，FBsnow 有解 %d。\n'], ...
                pairIndex, n, 100 * pairIndex / n, ...
                char(string(inputs.day(pairIndex), 'yyyy-MM-dd')), ...
                elapsedSeconds / 60, etaSeconds / 60, ...
                nnz(nRadarRoots(1:pairIndex) > 0), ...
                nnz(nLaserRoots(1:pairIndex) > 0));
        end
    end

    results = addvars(inputs, hiRadar, hsRadar, hiLaser, hsLaser, ...
        hiLaserAlt, hsLaserAlt, nRadarRoots, nLaserRoots, tbFromOIB);
    nearTB = abs(results.tbFromOIB - results.TB) <= 1.5;
    sortedCounts = sort(results.nOIB(isfinite(results.nOIB)));
    coverageThreshold = sortedCounts( ...
        max(1, ceil(0.95 * numel(sortedCounts))));
    highCoverage = results.nOIB >= coverageThreshold;
    metricSummary = summarize_retrieval_subset( ...
        results, true(n, 1), nearTB);
    highCoverageFit = fit_retrieval( ...
        results.TB(highCoverage), results.tbFromOIB(highCoverage));
    highCoverageRow = metricSummary(1, :);
    highCoverageRow.metric = "TB：OIB 样本数前 5%";
    fitFields = {'n','correlation','rSquared','bias','rmse','slope','intercept'};
    for fieldIndex = 1:numel(fitFields)
        fieldName = fitFields{fieldIndex};
        highCoverageRow.(fieldName) = highCoverageFit.(fieldName);
    end
    metricSummary = [metricSummary(1, :); highCoverageRow; ...
        metricSummary(2:end, :)];
    metricByYear = summarize_retrieval_by_year(results, nearTB, ...
        fullfile(resultRoot, outputStem + "_metrics_by_year.csv"));

    save(outputFile, 'pairs', 'summary', 'results', 'metricSummary', ...
        'metricByYear', 'nearTB', 'highCoverage', 'coverageThreshold', ...
        'model', 'operatorOptions', 'gridVersion', 'referenceMatch', ...
        'retrievalConfig', '-v7.3');
    writetable(results, fullfile(resultRoot, ...
        outputStem + "_retrieval.csv"));
    writetable(metricSummary, fullfile(resultRoot, ...
        outputStem + "_metrics.csv"));

    contextLabel = '步长 3，互不重叠 3×3 窗口，先平均 hi/hs 后计算 TB';
    plot_retrieval_validation(results, nearTB, highCoverage, ...
        figureRoot, contextLabel);
    fprintf('完成 %d 个飞行日、%d 条配对；逐条反演 %.1f min。\n', ...
        height(summary), n, toc(solverTimer) / 60);
    fprintf('结果：%s\n图件：%s\n', outputFile, figureRoot);
    disp(metricSummary);
end


function flightFiles = find_flight_files(dataRoot, targetDays)
    specifications = { ...
        fullfile(dataRoot, 'OIB', 'IDCSI4_L4', '2012'), ...
        'IDCSI4_*.txt'; ...
        fullfile(dataRoot, 'OIB', 'NSIDC-0708_QuickLook', ...
        '2015_GR_NASA'), 'OIB_*_IDCSI2*.txt'};
    days = datetime.empty(0, 1);
    paths = strings(0, 1);
    for specificationIndex = 1:size(specifications, 1)
        files = dir(fullfile(specifications{specificationIndex, 1}, ...
            specifications{specificationIndex, 2}));
        for fileIndex = 1:numel(files)
            token = regexp(files(fileIndex).name, ...
                '(\d{8})', 'match', 'once');
            day = datetime(token, 'InputFormat', 'yyyyMMdd');
            if isempty(targetDays) || ismember(day, targetDays)
                days(end + 1, 1) = day; %#ok<AGROW>
                paths(end + 1, 1) = string(fullfile( ...
                    files(fileIndex).folder, files(fileIndex).name)); %#ok<AGROW>
            end
        end
    end
    flightFiles = sortrows(table(days, paths, ...
        'VariableNames', {'day','path'}), 'day');
end


function referenceMatch = verify_against_stride_filter( ...
        pairs, matlabRoot, targetDays, operatorOptions)
%VERIFY_AGAINST_STRIDE_FILTER 确认本轮只改变中心格抽样步长。

    referenceMatch = false;
    referenceFile = fullfile(matlabRoot, 'processed', ...
        'pairs_2012_2015_same_day.mat');
    if ~isfile(referenceFile)
        warning('run_stride3_parameter_mean_inversion:NoReference', ...
            '缺少历史全中心格缓存，跳过逐值对照：%s', referenceFile);
        return;
    end
    referenceData = load(referenceFile, 'pairs');
    reference = referenceData.pairs;
    selected = mod(reference.row - operatorOptions.centerPhase(1), ...
        operatorOptions.centerStride) == 0 & ...
        mod(reference.col - operatorOptions.centerPhase(2), ...
        operatorOptions.centerStride) == 0;
    reference = reference(selected, :);
    if ~isempty(targetDays)
        reference = reference(ismember(reference.day, targetDays), :);
    end
    reference = sortrows(reference, {'day','row','col'});
    comparisonVariables = intersect(pairs.Properties.VariableNames, ...
        reference.Properties.VariableNames, 'stable');
    if ~isequaln(pairs(:, comparisonVariables), ...
            reference(:, comparisonVariables))
        error('run_stride3_parameter_mean_inversion:ReferenceMismatch', ...
            ['步长 3 重新配准结果与历史全中心格结果按相位筛选后不一致；' ...
            '不能把差异只归因于步长。']);
    end
    referenceMatch = true;
    fprintf('逐值核对通过：步长 3 配对等于历史配对的 [2,2] 相位子集。\n');
end
