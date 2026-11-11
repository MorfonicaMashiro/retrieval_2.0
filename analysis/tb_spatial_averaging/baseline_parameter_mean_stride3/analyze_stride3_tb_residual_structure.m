function report = analyze_stride3_tb_residual_structure()
%ANALYZE_STRIDE3_TB_RESIDUAL_STRUCTURE 拆解冻结步长 3 基线 TB 残差来源。
%   残差定义为 T_B,sim - T_B,SMOS。输出年份、冰型、覆盖度和飞行日
%   分组统计，列出 |残差|>=10 K 的长尾点，并绘制四面板诊断图。

    experimentRoot = fileparts(mfilename('fullpath'));
    matlabRoot = fileparts(fileparts(fileparts(experimentRoot)));
    resultRoot = fullfile(experimentRoot, 'results');
    figureRoot = fullfile(experimentRoot, 'figures');
    addpath(fullfile(matlabRoot, 'plot'));
    if ~isfolder(figureRoot), mkdir(figureRoot); end

    sourceFile = fullfile(resultRoot, ...
        'stride3_parameter_mean_2012_2015_retrieval.csv');
    paperFile = fullfile(matlabRoot, 'processed', ...
        'tb_forward_paper_baseline_v1_uncapped_table1.csv');
    if ~isfile(sourceFile)
        error('analyze_stride3_tb_residual_structure:MissingResult', ...
            '缺少冻结步长 3 基线结果：%s', sourceFile);
    end
    results = readtable(sourceFile);
    valid = isfinite(results.TB) & isfinite(results.tbFromOIB) & ...
        isfinite(results.nOIB);
    results = results(valid,:);
    results.differenceK = results.tbFromOIB - results.TB;
    results.absDifferenceK = abs(results.differenceK);
    results.observationYear = year(datetime(results.day));
    results.iceType = repmat("FYI", height(results), 1);
    results.iceType(logical(results.isMYI)) = "MYI";

    coverageEdges = prctile(results.nOIB, [25 50 75 95]);
    coverageGroup = strings(height(results), 1);
    coverageGroup(results.nOIB <= coverageEdges(1)) = "Q1最低25%";
    coverageGroup(results.nOIB > coverageEdges(1) & ...
        results.nOIB <= coverageEdges(2)) = "Q2";
    coverageGroup(results.nOIB > coverageEdges(2) & ...
        results.nOIB <= coverageEdges(3)) = "Q3";
    coverageGroup(results.nOIB > coverageEdges(3) & ...
        results.nOIB < coverageEdges(4)) = "Q4：75%–95%";
    coverageGroup(results.nOIB >= coverageEdges(4)) = "前5%";
    coverageOrder = ["Q1最低25%","Q2","Q3","Q4：75%–95%","前5%"];
    results.coverageGroup = categorical(coverageGroup, ...
        coverageOrder, 'Ordinal', true);

    groupSummary = summary_row("全部", "全部", results.differenceK);
    years = unique(results.observationYear).';
    for currentYear = years
        selected = results.observationYear == currentYear;
        groupSummary = [groupSummary; summary_row("年份", ...
            string(currentYear), results.differenceK(selected))]; %#ok<AGROW>
    end
    for iceType = ["FYI","MYI"]
        selected = results.iceType == iceType;
        groupSummary = [groupSummary; summary_row("冰型", ...
            iceType, results.differenceK(selected))]; %#ok<AGROW>
    end
    for currentYear = years
        for iceType = ["FYI","MYI"]
            selected = results.observationYear == currentYear & ...
                results.iceType == iceType;
            groupSummary = [groupSummary; summary_row("年份×冰型", ...
                string(currentYear) + " " + iceType, ...
                results.differenceK(selected))]; %#ok<AGROW>
        end
    end

    coverageSummary = table();
    for group = coverageOrder
        selected = string(results.coverageGroup) == group;
        current = summary_row("OIB覆盖组", group, ...
            results.differenceK(selected));
        current.minimumNOIB = min(results.nOIB(selected));
        current.maximumNOIB = max(results.nOIB(selected));
        current.medianNOIB = median(results.nOIB(selected));
        coverageSummary = [coverageSummary; current]; %#ok<AGROW>
    end

    uniqueDays = unique(results.day);
    dailySummary = table();
    for dayIndex = 1:numel(uniqueDays)
        currentDay = uniqueDays(dayIndex);
        selected = results.day == currentDay;
        current = summary_row("飞行日", string(currentDay), ...
            results.differenceK(selected));
        current.day = currentDay;
        dailySummary = [dailySummary; current]; %#ok<AGROW>
    end
    dailySummary = sortrows(dailySummary, 'meanDifferenceK', 'descend');

    tailThresholdK = 10;
    tailRows = results(results.absDifferenceK >= tailThresholdK, ...
        {'day','flight','row','col','TB','tbFromOIB','differenceK', ...
        'absDifferenceK','nOIB','isMYI','hi','hs'});
    tailRows = sortrows(tailRows, 'absDifferenceK', 'descend');

    q01 = prctile(results.differenceK, 1);
    q99 = prctile(results.differenceK, 99);
    robustSummary = [ ...
        summary_row("稳健性", "全部", results.differenceK); ...
        summary_row("稳健性", "去除|差值|≥10 K", ...
            results.differenceK(results.absDifferenceK < tailThresholdK)); ...
        summary_row("稳健性", "保留1%–99%", ...
            results.differenceK(results.differenceK >= q01 & ...
            results.differenceK <= q99))];

    paperSummary = table();
    if isfile(paperFile)
        paperCases = readtable(paperFile);
        paperSummary = summary_row("论文八场景", "全部", ...
            paperCases.differenceK);
        paperSummary.modelHigherCount = nnz(paperCases.differenceK > 0);
        paperSummary.paperHigherCount = nnz(paperCases.differenceK < 0);
        paperSummary.equalCount = nnz(paperCases.differenceK == 0);
    end

    [~, commonEdges] = histcounts(results.differenceK, 'BinMethod', 'fd');
    figureFile = fullfile(figureRoot, ...
        'tb_residual_structure_by_year_ice_coverage.png');
    draw_diagnostic_figure(results, commonEdges, figureFile);

    writetable(groupSummary, fullfile(resultRoot, ...
        'stride3_tb_residual_group_summary.csv'));
    writetable(coverageSummary, fullfile(resultRoot, ...
        'stride3_tb_residual_coverage_summary.csv'));
    writetable(dailySummary, fullfile(resultRoot, ...
        'stride3_tb_residual_daily_summary.csv'));
    writetable(tailRows, fullfile(resultRoot, ...
        'stride3_tb_residual_abs_ge_10K.csv'));
    writetable(robustSummary, fullfile(resultRoot, ...
        'stride3_tb_residual_robustness_summary.csv'));
    if ~isempty(paperSummary)
        writetable(paperSummary, fullfile(resultRoot, ...
            'table1_model_minus_paper_summary.csv'));
    end

    report = struct('groupSummary', groupSummary, ...
        'coverageSummary', coverageSummary, 'dailySummary', dailySummary, ...
        'tailRows', tailRows, 'robustSummary', robustSummary, ...
        'paperSummary', paperSummary, 'figureFile', string(figureFile), ...
        'coverageEdges', coverageEdges, 'tailThresholdK', tailThresholdK);
    fprintf('残差结构诊断图：%s\n', figureFile);
    disp(groupSummary);
    disp(coverageSummary);
    disp(robustSummary);
    if ~isempty(paperSummary), disp(paperSummary); end
end


function output = summary_row(groupType, groupLabel, values)
    values = double(values(:));
    values = values(isfinite(values));
    n = numel(values);
    if n == 0
        meanDifferenceK = NaN; medianDifferenceK = NaN;
        standardDeviationK = NaN; rmseK = NaN;
        percentile5K = NaN; percentile95K = NaN;
        positivePercent = NaN; within1p5KPercent = NaN;
    else
        meanDifferenceK = mean(values);
        medianDifferenceK = median(values);
        standardDeviationK = std(values);
        rmseK = sqrt(mean(values .^ 2));
        percentile5K = prctile(values, 5);
        percentile95K = prctile(values, 95);
        positivePercent = 100 * mean(values > 0);
        within1p5KPercent = 100 * mean(abs(values) <= 1.5);
    end
    output = table(string(groupType), string(groupLabel), n, ...
        meanDifferenceK, medianDifferenceK, standardDeviationK, rmseK, ...
        percentile5K, percentile95K, positivePercent, ...
        within1p5KPercent, 'VariableNames', {'groupType','group', ...
        'n','meanDifferenceK','medianDifferenceK','standardDeviationK', ...
        'rmseK','percentile5K','percentile95K','positivePercent', ...
        'within1p5KPercent'});
end


function draw_diagnostic_figure(results, edges, pngPath)
    colors = [0.00 0.45 0.70; 0.90 0.62 0.00; ...
        0.00 0.62 0.45; 0.80 0.20 0.15];
    centers = edges(1:end-1) + diff(edges) ./ 2;
    figureHandle = figure('Visible', 'off', 'Color', 'w', ...
        'Position', [50 50 1250 900]);
    layout = tiledlayout(figureHandle, 2, 2, ...
        'TileSpacing', 'compact', 'Padding', 'compact');
    title(layout, '步长 3 基线：模拟亮温减 SMOS 观测亮温的残差结构', ...
        'FontName', 'Microsoft YaHei', 'FontWeight', 'bold');

    ax1 = nexttile(layout);
    hold(ax1, 'on');
    years = unique(results.observationYear).';
    for index = 1:numel(years)
        selected = results.observationYear == years(index);
        frequency = 100 .* histcounts(results.differenceK(selected), edges) ...
            ./ nnz(selected);
        stairs(ax1, centers, frequency, 'LineWidth', 2.0, ...
            'Color', colors(index,:), ...
            'DisplayName', sprintf('%d（n=%d）', years(index), nnz(selected)));
    end
    style_histogram_axes(ax1, edges, '按年份');
    legend(ax1, 'Location', 'northeast');

    ax2 = nexttile(layout);
    hold(ax2, 'on');
    iceTypes = ["FYI","MYI"];
    for index = 1:numel(iceTypes)
        selected = results.iceType == iceTypes(index);
        frequency = 100 .* histcounts(results.differenceK(selected), edges) ...
            ./ nnz(selected);
        stairs(ax2, centers, frequency, 'LineWidth', 2.0, ...
            'Color', colors(index+2,:), ...
            'DisplayName', sprintf('%s（n=%d）', iceTypes(index), nnz(selected)));
    end
    style_histogram_axes(ax2, edges, '按冰型');
    legend(ax2, 'Location', 'northeast');

    ax3 = nexttile(layout);
    hold(ax3, 'on');
    fyi = results.iceType == "FYI";
    myi = results.iceType == "MYI";
    scatter(ax3, results.nOIB(fyi), results.differenceK(fyi), 18, ...
        colors(1,:), '^', 'filled', 'MarkerFaceAlpha', 0.42, ...
        'MarkerEdgeAlpha', 0.3, 'DisplayName', 'FYI');
    scatter(ax3, results.nOIB(myi), results.differenceK(myi), 18, ...
        colors(2,:), 's', 'filled', 'MarkerFaceAlpha', 0.42, ...
        'MarkerEdgeAlpha', 0.3, 'DisplayName', 'MYI');
    yline(ax3, 0, '-', 'Color', colors(4,:), 'LineWidth', 1.5, ...
        'HandleVisibility', 'off');
    set(ax3, 'XScale', 'log');
    xlabel(ax3, '单个配对的 OIB 样本数（对数轴）');
    ylabel(ax3, '差值（K）');
    title(ax3, '残差与 OIB 覆盖点数');
    legend(ax3, 'Location', 'southwest');
    style_axes(ax3);

    ax4 = nexttile(layout);
    boxchart(ax4, results.coverageGroup, results.differenceK, ...
        'BoxFaceColor', colors(1,:), 'BoxFaceAlpha', 0.62, ...
        'MarkerColor', colors(2,:));
    hold(ax4, 'on');
    yline(ax4, 0, '-', 'Color', colors(4,:), 'LineWidth', 1.5);
    xlabel(ax4, 'OIB 样本数分组');
    ylabel(ax4, '差值（K）');
    title(ax4, sprintf('按覆盖度分组（阈值：%g / %g / %g / %g 点）', ...
        prctile(results.nOIB,25), prctile(results.nOIB,50), ...
        prctile(results.nOIB,75), prctile(results.nOIB,95)));
    style_axes(ax4);
    ax4.XTickLabelRotation = 15;

    exportgraphics(figureHandle, pngPath, 'Resolution', 300);
    close(figureHandle);
end


function style_histogram_axes(ax, edges, panelTitle)
    colors = [0.80 0.20 0.15];
    xline(ax, 0, '-', 'Color', colors, 'LineWidth', 1.4, ...
        'HandleVisibility', 'off');
    xlim(ax, [edges(1), edges(end)]);
    xlabel(ax, '差值 T_{B,sim} - T_{B,SMOS}（K）');
    ylabel(ax, '组内频率（%）');
    title(ax, panelTitle);
    style_axes(ax);
end


function style_axes(ax)
    grid(ax, 'on');
    box(ax, 'on');
    ax.GridAlpha = 0.18;
    ax.FontName = 'Microsoft YaHei';
    ax.FontSize = 10;
end
