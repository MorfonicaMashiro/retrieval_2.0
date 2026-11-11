function [pairs, summary, metrics, differences] = ...
        run_sample_tb_mean_operator(targetDays, outputStem)
%RUN_SAMPLE_TB_MEAN_OPERATOR 运行逐 OIB 点算 TB 后平均的前向试验。
%   TARGETDAYS 为空时处理 2012、2015 年全部本地 OIB 飞行日；也可以传入
%   datetime 向量只处理指定日期。OUTPUTSTEM 是 results/ 下的输出前缀。
%   表温和 Vb 截断口径由 retrieval_config 统一控制。

    if nargin < 1
        targetDays = datetime.empty(0, 1);
    end
    targetDays = dateshift(datetime(targetDays(:)), 'start', 'day');

    operatorRoot = fileparts(mfilename('fullpath'));
    matlabRoot = fileparts(fileparts(fileparts(operatorRoot)));
    dataRoot = fullfile(fileparts(matlabRoot), 'data');
    resultRoot = fullfile(operatorRoot, 'results');
    addpath(matlabRoot, fullfile(matlabRoot, 'readdata'), operatorRoot);
    if ~isfolder(resultRoot), mkdir(resultRoot); end

    retrievalConfig = retrieval_config();
    retrievalConfig.collocation.computeSampleMeanTB = true;
    retrievalConfig.collocation.parameterWeighting = "none";
    model = retrievalConfig.model;
    gridVersion = retrievalConfig.data.gridVersion;
    if nargin < 2 || strlength(string(outputStem)) == 0
        outputStem = "sample_tb_mean_" + ...
            retrievalConfig.output.variantTag + "_2012_2015";
    end
    outputStem = string(outputStem);
    flightFiles = find_flight_files(dataRoot, targetDays);
    if isempty(flightFiles)
        error('run_sample_tb_mean_operator:NoFlights', ...
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
            error('run_sample_tb_mean_operator:MissingSMOS', ...
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
            stats.mixedTypeCells, stats.forwardPoints, ...
            stats.forwardSeconds, ...
            'VariableNames', {'day','rawPoints','qualityPoints', ...
            'tbOverlapCells','matchedCells','mixedTypeCells', ...
            'forwardPoints','forwardSeconds'})]; %#ok<AGROW>

        save(outputFile, 'pairs', 'summary', 'model', ...
            'gridVersion', 'retrievalConfig', '-v7.3');
        elapsedSeconds = toc(timer);
        etaSeconds = elapsedSeconds / flightIndex * ...
            (height(flightFiles) - flightIndex);
        fprintf(['[%d/%d, %5.1f%%] %s：逐点前向 %d 个，配对 %d 条；' ...
            '本日前向 %.1f min，累计 %.1f min，预计剩余 %.1f min。\n'], ...
            flightIndex, height(flightFiles), ...
            100 * flightIndex / height(flightFiles), token, ...
            stats.forwardPoints, stats.matchedCells, ...
            stats.forwardSeconds / 60, elapsedSeconds / 60, ...
            etaSeconds / 60);
    end

    [~, order] = sortrows([datenum(pairs.day), pairs.row, pairs.col]); %#ok<DATNM>
    pairs = pairs(order, :);
    summary = sortrows(summary, 'day');
    [metrics, differences] = summarize_sample_tb_mean_operator(pairs);
    save(outputFile, 'pairs', 'summary', 'metrics', 'differences', ...
        'model', 'gridVersion', 'retrievalConfig', '-v7.3');
    writetable(pairs, fullfile(resultRoot, outputStem + "_pairs.csv"));
    writetable(summary, fullfile(resultRoot, outputStem + "_summary.csv"));
    writetable(metrics, fullfile(resultRoot, outputStem + "_metrics.csv"));
    writetable(differences, ...
        fullfile(resultRoot, outputStem + "_differences.csv"));
    fprintf('完成 %d 个飞行日、%d 条配对：%s\n', ...
        height(summary), height(pairs), outputFile);
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
            token = regexp(files(fileIndex).name, '(\d{8})', 'match', 'once');
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
