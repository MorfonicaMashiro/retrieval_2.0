function [results, metricSummary, metricByYear] = ...
        run_kt19_ease2_inversion(forwardStem, outputStem)
%RUN_KT19_EASE2_INVERSION 反演 KT19--EASE2 配对的冰厚与雪深。
%   对每个配对分别使用其窗口平均 KT19_surf，运行 TB+FBice 和
%   TB+FBsnow 两条路线。输入、结果和图件均留在本试验分支，不覆盖
%   历史 EASE 结果。

    if nargin < 1 || strlength(string(forwardStem)) == 0
        forwardStem = "kt19_ease2_linear_extrapolation_2012_2015";
    end
    if nargin < 2 || strlength(string(outputStem)) == 0
        outputStem = "kt19_ease2_linear_extrapolation_2012_2015_inversion";
    end
    forwardStem = string(forwardStem);
    outputStem = string(outputStem);

    branchRoot = fileparts(mfilename('fullpath'));
    matlabRoot = fileparts(fileparts(fileparts(branchRoot)));
    resultRoot = fullfile(branchRoot, 'results');
    figureRoot = fullfile(branchRoot, 'figures', 'inversion');
    addpath(matlabRoot, fullfile(matlabRoot, 'plot'), branchRoot);
    if ~isfolder(resultRoot), mkdir(resultRoot); end
    if ~isfolder(figureRoot), mkdir(figureRoot); end

    forwardFile = fullfile(resultRoot, forwardStem + ".mat");
    if ~isfile(forwardFile)
        error('run_kt19_ease2_inversion:MissingForwardResult', ...
            '缺少前向配对结果：%s', forwardFile);
    end
    source = load(forwardFile, 'pairs', 'summary', 'config', 'model');
    pairs = source.pairs;
    summary = source.summary;
    config = source.config;
    model = source.model;
    if string(config.data.gridVersion) ~= "EASE2" || ...
            string(config.temperature.source) ~= "oib_kt19" || ...
            ~ismember('surfaceTempC', pairs.Properties.VariableNames)
        error('run_kt19_ease2_inversion:Configuration', ...
            '输入必须是含逐配对 KT19 表温的 EASE2 结果。');
    end

    inputNames = {'day','flight','row','col','TB','FBice', ...
        'FBsnow','isMYI','hi','hs','nOIB','surfaceTempC', ...
        'TBOIBKT19Linear'};
    inputs = pairs(:, inputNames);
    inputs.Properties.VariableNames{'TBOIBKT19Linear'} = 'tbFromOIB';
    n = height(inputs);
    hiRadar = nan(n,1); hsRadar = nan(n,1);
    hiLaser = nan(n,1); hsLaser = nan(n,1);
    hiLaserAlt = nan(n,1); hsLaserAlt = nan(n,1);
    nRadarRoots = zeros(n,1); nLaserRoots = zeros(n,1);
    radarStatus = strings(n,1); laserStatus = strings(n,1);
    outputFile = fullfile(resultRoot, outputStem + ".mat");
    solverTimer = tic;
    progressEvery = max(1,ceil(n/50));
    fprintf(['开始 KT19--EASE2 逐对反演：%d 对；' ...
        'TB+FBice 与 TB+FBsnow。\n'],n);
    for pairIndex = 1:n
        pairConfig = config;
        pairConfig.model = model;
        pairConfig.model.surfaceTempC = inputs.surfaceTempC(pairIndex);
        [radar, radarInfo] = invert_pair( ...
            inputs.TB(pairIndex), inputs.FBice(pairIndex), "ice", ...
            inputs.isMYI(pairIndex), pairConfig);
        [laser, laserInfo] = invert_pair( ...
            inputs.TB(pairIndex), inputs.FBsnow(pairIndex), "snow", ...
            inputs.isMYI(pairIndex), pairConfig);
        radarStatus(pairIndex) = radarInfo.status;
        laserStatus(pairIndex) = laserInfo.status;
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

        if pairIndex == 1 || mod(pairIndex,progressEvery) == 0 || ...
                pairIndex == n
            elapsedMinutes = toc(solverTimer)/60;
            remainingMinutes = elapsedMinutes/pairIndex*(n-pairIndex);
            fprintf(['[%d/%d, %5.1f%%] 有解：FBice %d，FBsnow %d；' ...
                '已用 %.1f min，预计剩余 %.1f min。\n'], ...
                pairIndex,n,100*pairIndex/n, ...
                nnz(nRadarRoots(1:pairIndex)>0), ...
                nnz(nLaserRoots(1:pairIndex)>0), ...
                elapsedMinutes,remainingMinutes);
        end
        if mod(pairIndex,50) == 0
            lastCompletedIndex = pairIndex;
            save(outputFile,'inputs','hiRadar','hsRadar','hiLaser', ...
                'hsLaser','hiLaserAlt','hsLaserAlt','nRadarRoots', ...
                'nLaserRoots','radarStatus','laserStatus', ...
                'lastCompletedIndex','config','model', ...
                '-v7.3');
        end
    end

    results = addvars(inputs,hiRadar,hsRadar,hiLaser,hsLaser, ...
        hiLaserAlt,hsLaserAlt,nRadarRoots,nLaserRoots, ...
        radarStatus,laserStatus);
    nearThresholdK = config.validation.nearTbThresholdK;
    nearTB = abs(results.tbFromOIB-results.TB) <= nearThresholdK;
    sortedCounts = sort(results.nOIB(isfinite(results.nOIB)));
    quantileIndex = max(1,ceil(config.validation.highCoverageQuantile * ...
        numel(sortedCounts)));
    coverageThreshold = sortedCounts(quantileIndex);
    highCoverage = results.nOIB >= coverageThreshold;
    metricSummary = summarize_retrieval_subset( ...
        results,true(n,1),nearTB,nearThresholdK);
    highCoverageFit = fit_retrieval( ...
        results.TB(highCoverage),results.tbFromOIB(highCoverage));
    highCoverageRow = metricSummary(1,:);
    highCoverageRow.metric = "TB：OIB 样本数前 5%";
    fitFields = {'n','correlation','rSquared','bias','rmse', ...
        'slope','intercept'};
    for fieldIndex = 1:numel(fitFields)
        fieldName = fitFields{fieldIndex};
        highCoverageRow.(fieldName) = highCoverageFit.(fieldName);
    end
    metricSummary = [metricSummary(1,:); highCoverageRow; ...
        metricSummary(2:end,:)];
    metricByYear = summarize_retrieval_by_year(results,nearTB, ...
        fullfile(resultRoot,outputStem + "_metrics_by_year.csv"), ...
        nearThresholdK);

    method = struct('grid',"EASE2-N 12.5 km", ...
        'temperature',"pair-specific mean OIB KT19_surf", ...
        'coldBrinePolicy', ...
        "Cox-Weeks F1 tangent-linear extrapolation below -30 C", ...
        'radiativeTransfer',"Burke", ...
        'routes',["TB+FBice","TB+FBsnow"]);
    save(outputFile,'pairs','summary','results','metricSummary', ...
        'metricByYear','nearTB','highCoverage','coverageThreshold', ...
        'config','model','method','-v7.3');
    writetable(results,fullfile(resultRoot,outputStem + "_retrieval.csv"));
    writetable(metricSummary,fullfile(resultRoot,outputStem + ...
        "_metrics.csv"));

    contextLabel = ['EASE2；步长3；逐配对 OIB KT19；' ...
        'Cox–Weeks 冷端线性外推'];
    plot_retrieval_validation(results,nearTB,highCoverage, ...
        figureRoot,contextLabel,nearThresholdK, ...
        config.validation.highCoverageQuantile);
    fprintf(['KT19--EASE2 反演完成：%d 对；FBice 有解 %d，' ...
        'FBsnow 有解 %d；%.1f min。\n%s\n'],n, ...
        nnz(nRadarRoots>0),nnz(nLaserRoots>0), ...
        toc(solverTimer)/60,outputFile);
    disp(metricSummary);
end
