function [summary, bins] = plot_stride3_tb_residual_distribution()
%PLOT_STRIDE3_TB_RESIDUAL_DISTRIBUTION 绘制冻结步长 3 基线 TB 残差分布。

    experimentRoot = fileparts(mfilename('fullpath'));
    matlabRoot = fileparts(fileparts(fileparts(experimentRoot)));
    resultRoot = fullfile(experimentRoot, 'results');
    figureRoot = fullfile(experimentRoot, 'figures');
    addpath(fullfile(matlabRoot, 'plot'));

    sourceFile = fullfile(resultRoot, ...
        'stride3_parameter_mean_2012_2015_retrieval.csv');
    if ~isfile(sourceFile)
        error('plot_stride3_tb_residual_distribution:MissingResult', ...
            '缺少冻结步长 3 基线结果：%s', sourceFile);
    end
    results = readtable(sourceFile);
    contextLabel = '2012 与 2015 年；步长 3；先平均 h_i/h_s 后计算亮温';
    figureFile = fullfile(figureRoot, ...
        'tb_simulated_minus_smos_frequency_distribution.png');
    [summary, bins] = plot_tb_residual_distribution( ...
        results, figureFile, contextLabel);
    writetable(summary, fullfile(resultRoot, ...
        'stride3_parameter_mean_2012_2015_tb_residual_summary.csv'));
    writetable(bins, fullfile(resultRoot, ...
        'stride3_parameter_mean_2012_2015_tb_residual_bins.csv'));
    fprintf('亮温差值分布图：%s\n', figureFile);
    disp(summary);
end
