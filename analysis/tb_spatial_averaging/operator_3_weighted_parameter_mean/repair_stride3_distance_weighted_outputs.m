function metricSummary = repair_stride3_distance_weighted_outputs(outputStem)
%REPAIR_STRIDE3_DISTANCE_WEIGHTED_OUTPUTS 修复前 5% 指标和对应 TB 图。
%   复用已完成的反演结果，不重新求解冰厚和雪深。

    if nargin < 1 || strlength(string(outputStem)) == 0
        outputStem = "stride3_distance_weighted_parameters_2012_2015";
    end
    outputStem = string(outputStem);
    routeRoot = fileparts(mfilename('fullpath'));
    matlabRoot = fileparts(fileparts(fileparts(routeRoot)));
    resultRoot = fullfile(routeRoot, 'results', 'inversion');
    figureRoot = fullfile(routeRoot, 'figures', 'inversion');
    addpath(matlabRoot, fullfile(matlabRoot, 'plot'));
    outputFile = fullfile(resultRoot, outputStem + ".mat");
    if ~isfile(outputFile)
        error('repair_stride3_distance_weighted_outputs:MissingResult', ...
            '没有找到全量反演结果：%s', outputFile);
    end

    saved = load(outputFile);
    pairs = saved.pairs;
    summary = saved.summary;
    results = saved.results;
    model = saved.model;
    operatorOptions = saved.operatorOptions;
    metricByYear = saved.metricByYear;
    if height(results) ~= height(pairs)
        error('repair_stride3_distance_weighted_outputs:PairCount', ...
            'results 与 pairs 的行数不一致。');
    end

    % 论文的“样本数前 5%”按原始 OIB 点数定义，不按权重支持数定义。
    results.nOIB = pairs.nOIB;
    if ismember('nOIBWeightedSupport', results.Properties.VariableNames)
        results.nOIBWeightedSupport = pairs.nOIBWeightedSupport;
    else
        results = addvars(results, pairs.nOIBWeightedSupport, ...
            'After', 'nOIB', 'NewVariableNames', 'nOIBWeightedSupport');
    end
    nearTB = abs(results.tbFromOIB - results.TB) <= 1.5;
    sortedCounts = sort(results.nOIB(isfinite(results.nOIB)));
    coverageThreshold = sortedCounts( ...
        max(1, ceil(0.95 * numel(sortedCounts))));
    highCoverage = results.nOIB >= coverageThreshold;

    metricSummary = summarize_retrieval_subset( ...
        results, true(height(results), 1), nearTB);
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

    save(outputFile, 'pairs', 'summary', 'results', 'metricSummary', ...
        'metricByYear', 'nearTB', 'highCoverage', 'coverageThreshold', ...
        'model', 'operatorOptions', '-v7.3');
    writetable(results, fullfile(resultRoot, ...
        outputStem + "_retrieval.csv"));
    writetable(metricSummary, fullfile(resultRoot, ...
        outputStem + "_metrics.csv"));
    plot_retrieval_validation(results, nearTB, highCoverage, figureRoot, ...
        '到 SMOS 中心反距离加权参数→单一冰雪柱→单次 TB');
    fprintf('前 5%% 阈值 nOIB >= %g，共 %d 条配对。\n', ...
        coverageThreshold, nnz(highCoverage));
end
