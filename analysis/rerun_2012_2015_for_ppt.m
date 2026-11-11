% 仅保留 2012 和 2015 年样本，复算论文可比口径的验证指标与图件。
% 逐样本反演彼此独立，因此从四年逐条求解结果中筛选这两年，
% 与只将这两年样本送入相同求解器得到的验证结果等价。

analysisDir = fileparts(mfilename('fullpath'));
matlabRoot = fileparts(analysisDir);
addpath(matlabRoot);
addpath(fullfile(matlabRoot, 'plot'));

sourceCsv = fullfile(matlabRoot, 'processed', 'inversion_2012_2015.csv');
if ~isfile(sourceCsv)
    error('rerun_2012_2015_for_ppt:MissingInput', ...
        '缺少四年逐条反演结果：%s', sourceCsv);
end

allResults = readtable(sourceCsv, 'TextType', 'string');
if ~isdatetime(allResults.day)
    allResults.day = datetime(allResults.day, 'InputFormat', 'yyyy-MM-dd');
end
targetYears = [2012 2015];
yearMask = ismember(year(allResults.day), targetYears);
results = allResults(yearMask, :);
if height(results) ~= 6542
    warning('rerun_2012_2015_for_ppt:UnexpectedCount', ...
        '2012 与 2015 年样本数为 %d（预期 6542）。', height(results));
end

nearTB = abs(results.tbFromOIB - results.TB) <= 1.5;
sortedCounts = sort(results.nOIB(isfinite(results.nOIB)));
coverageThreshold = sortedCounts(max(1, ceil(0.95 * numel(sortedCounts))));
highCoverage = results.nOIB >= coverageThreshold;

metrics = summarize_retrieval_subset( ...
    results, true(height(results), 1), nearTB);
outputDir = fullfile(matlabRoot, 'figures', 'inversion_2012_2015_only');
if ~isfolder(outputDir), mkdir(outputDir); end
metricsCsv = fullfile(outputDir, 'metrics_2012_2015_only.csv');
resultsCsv = fullfile(outputDir, 'results_2012_2015_only.csv');
writetable(metrics, metricsCsv);
writetable(results, resultsCsv);
save(fullfile(outputDir, 'results_2012_2015_only.mat'), ...
    'results', 'metrics', 'nearTB', 'highCoverage', ...
    'coverageThreshold', 'targetYears', '-v7.3');

plot_retrieval_validation(results, nearTB, highCoverage, outputDir);

fprintf('仅 2012 与 2015 年：%d 条同日配对；1.5 K 筛选保留 %d 条。\n', ...
    height(results), nnz(nearTB));
disp(metrics);
fprintf('图件与指标已写入：%s\n', outputDir);
