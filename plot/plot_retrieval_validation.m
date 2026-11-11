function plot_retrieval_validation( ...
        results, nearTB, highCoverage, outputDir, contextLabel, ...
        nearTbThresholdK, highCoverageQuantile)
%PLOT_RETRIEVAL_VALIDATION 绘制论文图 11/12 风格的反演验证与亮温比较。
%   横轴为 OIB/SMOS 观测，纵轴为反演量/模型亮温；FYI 与 MYI 用不同符号。
%   CONTEXTLABEL 为可选实验说明；为空时保持历史图题不变。

    arguments
        results table
        nearTB (:,1) logical
        highCoverage (:,1) logical
        outputDir {mustBeTextScalar}
        contextLabel {mustBeTextScalar} = ""
        nearTbThresholdK (1,1) double {mustBePositive} = 1.5
        highCoverageQuantile (1,1) double {mustBeGreaterThan(highCoverageQuantile,0), ...
            mustBeLessThan(highCoverageQuantile,1)} = 0.95
    end
    if height(results) ~= numel(nearTB) || height(results) ~= numel(highCoverage)
        error('plot_retrieval_validation:MaskSize', '筛选掩膜长度不一致。');
    end
    outputDir = char(string(outputDir));
    contextLabel = char(string(contextLabel));
    if ~isfolder(outputDir), mkdir(outputDir); end

    resultYears = unique(year(results.day)).';
    if isequal(resultYears, [2012 2015])
        allPairsHeading = '2012 与 2015 年全部有效配对';
    elseif numel(resultYears) == 1
        allPairsHeading = sprintf('%d 年全部有效配对', resultYears);
    else
        allPairsHeading = sprintf('%d–%d 年全部有效配对', ...
            min(resultYears), max(resultYears));
    end
    if ~isempty(contextLabel)
        allPairsHeading = sprintf('%s\n%s', ...
            allPairsHeading, contextLabel);
    end
    plot_four_panel(results, true(height(results),1), ...
        fullfile(outputDir, 'retrieval_validation_all_revised.png'), ...
        allPairsHeading);
    thresholdTag = strrep(sprintf('%g', nearTbThresholdK), '.', 'p');
    plot_four_panel(results, nearTB, fullfile(outputDir, ...
        ['retrieval_validation_tb_within_' thresholdTag 'K_revised.png']), ...
        make_heading(sprintf( ...
        '|T_{B,model}-T_{B,SMOS}| ≤ %g K', nearTbThresholdK), ...
        contextLabel));
    plot_tb_comparison(results, highCoverage, ...
        fullfile(outputDir, 'tb_model_vs_smos_revised.png'), ...
        contextLabel, highCoverageQuantile);
end


function heading = make_heading(baseHeading, contextLabel)
    heading = baseHeading;
    if ~isempty(contextLabel)
        heading = sprintf('%s\n%s', baseHeading, contextLabel);
    end
end


function plot_four_panel(results, mask, pngPath, heading)
    observed = {results.hi, results.hs, results.hi, results.hs};
    retrieved = {results.hiRadar, results.hsRadar, ...
        results.hiLaser, results.hsLaser};
    titles = {'TB + FBice：冰厚', 'TB + FBice：雪深', ...
        'TB + FBsnow：冰厚', 'TB + FBsnow：雪深'};
    labels = {'冰厚 (m)', '雪深 (m)', '冰厚 (m)', '雪深 (m)'};
    figureHandle = figure('Visible', 'off', 'Color', 'w', ...
        'Position', [80 80 1120 900]);
    layout = tiledlayout(figureHandle, 2, 2, ...
        'TileSpacing', 'compact', 'Padding', 'compact');
    title(layout, heading, 'FontWeight', 'bold');
    for panel = 1:4
        ax = nexttile(layout);
        x = observed{panel}; y = retrieved{panel};
        valid = mask & isfinite(x) & isfinite(y);
        draw_scatter(ax, x, y, results.isMYI, valid, labels{panel}, ...
            titles{panel}, true);
    end
    export_figure(figureHandle, pngPath);
    close(figureHandle);
end


function plot_tb_comparison( ...
        results, highCoverage, pngPath, contextLabel, highCoverageQuantile)
    figureHandle = figure('Visible', 'off', 'Color', 'w', ...
        'Position', [80 80 1120 500]);
    layout = tiledlayout(figureHandle, 1, 2, ...
        'TileSpacing', 'compact', 'Padding', 'compact');
    if ~isempty(contextLabel)
        title(layout, contextLabel, 'FontWeight', 'bold');
    end
    masks = {true(height(results),1), highCoverage};
    highCoveragePercent = 100 * (1 - highCoverageQuantile);
    titles = {'全部有效配对', ...
        sprintf('OIB 样本数前 %g%%', highCoveragePercent)};
    for panel = 1:2
        ax = nexttile(layout);
        valid = masks{panel} & isfinite(results.TB) & isfinite(results.tbFromOIB);
        draw_scatter(ax, results.TB, results.tbFromOIB, ...
            results.isMYI, valid, '亮温 (K)', titles{panel}, false);
        xlabel(ax, 'SMOS 观测亮温 (K)');
        ylabel(ax, '由 OIB 冰厚/雪深模拟的亮温 (K)');
    end
    export_figure(figureHandle, pngPath);
    close(figureHandle);
end


function export_figure(figureHandle, pngPath)
%EXPORT_FIGURE 输出适合论文与汇报使用的 300 dpi PNG。
    exportgraphics(figureHandle, pngPath, 'Resolution', 300);
end


function draw_scatter(ax, x, y, isMYI, valid, axisLabel, panelTitle, includeZero)
    hold(ax, 'on');
    fyi = valid & ~isMYI;
    myi = valid & isMYI;
    scatter(ax, x(fyi), y(fyi), 17, [0.00 0.45 0.70], '^', ...
        'filled', 'MarkerFaceAlpha', 0.45, 'MarkerEdgeAlpha', 0.35, ...
        'DisplayName', 'FYI');
    scatter(ax, x(myi), y(myi), 17, [0.90 0.62 0.00], 's', ...
        'filled', 'MarkerFaceAlpha', 0.45, 'MarkerEdgeAlpha', 0.35, ...
        'DisplayName', 'MYI');
    finiteValues = [x(valid); y(valid)];
    if isempty(finiteValues)
        limits = [0 1];
    else
        lower = min(finiteValues);
        if includeZero, lower = min(0, lower); end
        upper = max(finiteValues);
        padding = max(0.02 * (upper - lower), eps);
        limits = [lower - padding, upper + padding];
    end
    plot(ax, limits, limits, 'k-', 'LineWidth', 1.0, ...
        'DisplayName', '1:1');
    fit = fit_retrieval(x(valid), y(valid));
    if isfinite(fit.slope)
        plot(ax, limits, fit.slope .* limits + fit.intercept, '--', ...
            'Color', [0.35 0.35 0.35], 'LineWidth', 1.1, ...
            'DisplayName', 'OLS');
    end
    if isfinite(fit.bias)
        plot(ax, limits, limits + fit.bias, ':', ...
            'Color', [0.75 0.20 0.15], 'LineWidth', 1.1, ...
            'DisplayName', '斜率=1');
    end
    xlim(ax, limits); ylim(ax, limits); axis(ax, 'square');
    grid(ax, 'on'); box(ax, 'on');
    xlabel(ax, ['OIB 观测' axisLabel]);
    ylabel(ax, ['反演' axisLabel]);
    title(ax, panelTitle);
    text(ax, 0.03, 0.97, sprintf(['n=%d\nR=%.2f, R^2=%.2f\n' ...
        'bias=%.3g, RMSE=%.3g'], fit.n, fit.correlation, ...
        fit.rSquared, fit.bias, fit.rmse), ...
        'Units', 'normalized', 'VerticalAlignment', 'top', ...
        'BackgroundColor', 'w', 'Margin', 3, 'FontSize', 9);
    if fit.n > 0, legend(ax, 'Location', 'southeast'); end
end
