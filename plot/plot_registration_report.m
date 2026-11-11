function outputs = plot_registration_report(outputDirectory)
%PLOT_REGISTRATION_REPORT 生成组会汇报使用的网格和配准质检图。
%   本函数不调用亮温前向模型，只使用原始 SMOS/OIB、EASE 网格定义和
%   已生成的同日配准缓存。输出图用于说明：
%     1) 同一天亮温在 EPSG:3413 源网格和原版北极 EASE 上的差异；
%     2) EASE 网格上的亮温、OIB 航迹和最终配对位置；
%     3) 23 个飞行日逐日匹配数量；
%     4) 6542 条配对样本的观测量和样本数分布。

    projectRoot = fileparts(fileparts(mfilename('fullpath')));
    dataRoot = fullfile(fileparts(projectRoot), 'data');
    addpath(projectRoot, fullfile(projectRoot, 'readdata'));
    if nargin < 1 || strlength(string(outputDirectory)) == 0
        outputDirectory = fullfile(projectRoot, 'figures', 'group_report');
    end
    outputDirectory = char(string(outputDirectory));
    if ~isfolder(outputDirectory), mkdir(outputDirectory); end

    outputs = struct();
    outputs.gridComparison = fullfile(outputDirectory, ...
        'grid_rotation_comparison.png');
    % 保留旧字段，避免已有脚本读取 outputs.gridRotation 时失效。
    outputs.gridRotation = outputs.gridComparison;
    outputs.registrationDay = fullfile(outputDirectory, ...
        'registration_20120314.png');
    outputs.dailySummary = fullfile(outputDirectory, ...
        'daily_matching_summary.png');
    outputs.pairOverview = fullfile(outputDirectory, ...
        'pair_sample_overview.png');

    fontName = 'Microsoft YaHei';
    green = [36 91 27] / 255;
    lightGreen = [151 191 125] / 255;
    grey = [0.72 0.74 0.75];
    orange = [217 120 33] / 255;
    blue = [49 112 167] / 255;

    plot_grid_tb_comparison(dataRoot, ...
        outputs.gridComparison, fontName, green);

    cachePath = fullfile(projectRoot, 'processed', ...
        'pairs_2017_same_day_v2.mat');
    cache = load(cachePath, 'pairs', 'summary');
    plot_registration_day(dataRoot, outputs.registrationDay, ...
        fontName, green);
    plot_daily_summary(cache.summary, outputs.dailySummary, ...
        fontName, green, grey, orange);
    plot_pair_overview(cache.pairs, outputs.pairOverview, ...
        fontName, green, lightGreen, orange, blue);
end

function plot_grid_tb_comparison(dataRoot, outputPath, ...
        fontName, green)
    token = '20120314';
    smosFile = fullfile(dataRoot, 'smos', '2012', ...
        ['SMOS_TB_north_' token '.nc']);
    smos = read_smos(smosFile);
    smosEase = grid_resampling.regrid_smos_to_ease(smos, "EASE");

    sourceCRS = projcrs(3413);
    [sourceX2D, sourceY2D] = projfwd(sourceCRS, ...
        double(smos.latitude), double(smos.longitude));
    % 原生 PSP 数据是 896×608（行×列）的矩形栅格。用每列/每行的
    % 投影坐标中位数恢复两条坐标轴，并保留完整矩阵，不能裁成方形。
    sourceX = median(sourceX2D, 1, 'omitnan') / 1000;
    sourceY = median(sourceY2D, 2, 'omitnan') / 1000;
    sourceTB = smos.tb;
    sourceTB(~isfinite(smos.nPair) | smos.nPair <= 0) = NaN;

    gridDef = smosEase.grid;
    easeX = gridDef.xCenters / 1000;
    easeY = gridDef.yCenters / 1000;
    fig = figure('Visible', 'off', 'Color', 'w', ...
        'Position', [80 80 1540 720]);
    tiles = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', ...
        'Padding', 'compact');

    ax1 = nexttile(tiles);
    h1 = imagesc(ax1, sourceX, sourceY, sourceTB, [205 255]);
    h1.AlphaData = isfinite(sourceTB);
    set(ax1, 'YDir', 'normal', 'Color', [0.94 0.95 0.94]);
    axis(ax1, 'equal', 'tight');
    xlabel(ax1, 'EPSG:3413 x (km)'); ylabel(ax1, 'EPSG:3413 y (km)');
    title(ax1, '原生 PSP：完整 896 × 608（行 × 列）', ...
        'FontWeight', 'bold');

    ax2 = nexttile(tiles);
    h2 = imagesc(ax2, easeX, easeY, smosEase.tb, [205 255]);
    h2.AlphaData = isfinite(smosEase.tb);
    set(ax2, 'YDir', 'normal', 'Color', [0.94 0.95 0.94]);
    axis(ax2, 'equal', 'tight');
    xlabel(ax2, 'EASE x (km)'); ylabel(ax2, 'EASE y (km)');
    title(ax2, '完整 12.5 km EASE：1441 × 1441', ...
        'FontWeight', 'bold');

    colormap(fig, parula(256));
    cb = colorbar(ax2, 'eastoutside');
    cb.Label.String = 'SMOS TB (K)';
    set([ax1 ax2], 'FontName', fontName, 'FontSize', 11.5, ...
        'Layer', 'top', 'Box', 'on');
    title(tiles, '2012-03-14 同日亮温在两套完整网格上的覆盖对比', ...
        'FontName', fontName, 'FontSize', 18, 'FontWeight', 'bold', ...
        'Color', green);
    exportgraphics(fig, outputPath, 'Resolution', 190);
    close(fig);
end

function plot_registration_day(dataRoot, outputPath, fontName, green)
    token = '20120314';
    smosFile = fullfile(dataRoot, 'smos', '2012', ...
        ['SMOS_TB_north_' token '.nc']);
    oibFile = fullfile(dataRoot, 'OIB', 'IDCSI4_L4', '2012', ...
        ['IDCSI4_' token '.txt']);
    smos = read_smos(smosFile);
    oib = read_oib(oibFile);
    smosEase = grid_resampling.regrid_smos_to_ease(smos, "EASE");
    [pairs, ~] = collocate_day(oib, smosEase, datetime(token, ...
        'InputFormat', 'yyyyMMdd'));

    tbl = oib.table;
    fbIce = double(tbl.mean_fb) - double(tbl.snow_depth);
    valid = isfinite(tbl.lat) & isfinite(tbl.lon) & ...
        isfinite(tbl.mean_fb) & isfinite(tbl.snow_depth) & ...
        isfinite(tbl.thickness) & isfinite(fbIce) & ...
        tbl.pcnt_ow == 0 & tbl.pcnt_thin_ice == 0 & ...
        tbl.pcnt_grey_ice == 0 & ismember(tbl.my_ice_flag, [0 1]);
    [oibRow, oibCol, oibX, oibY] = grid_resampling.PSP2EASE( ...
        tbl.lat(valid), tbl.lon(valid), "EASE");
    located = isfinite(oibRow) & isfinite(oibCol);
    oibX = oibX(located) / 1000;
    oibY = oibY(located) / 1000;

    gridDef = smosEase.grid;
    pairX = gridDef.xCenters(pairs.col) / 1000;
    pairY = gridDef.yCenters(pairs.row) / 1000;

    % 航迹总体近南北向。采用方形显示范围，避免等比例坐标下图面过窄。
    trackX = [oibX(:); pairX(:)];
    trackY = [oibY(:); pairY(:)];
    span = max(max(trackX) - min(trackX), ...
        max(trackY) - min(trackY)) + 350;
    span = max(span, 2200);
    centerX = mean([min(trackX), max(trackX)]);
    centerY = mean([min(trackY), max(trackY)]);
    xLimits = centerX + [-0.5 0.5] * span;
    yLimits = centerY + [-0.5 0.5] * span;
    xCentersKm = gridDef.xCenters / 1000;
    yCentersKm = gridDef.yCenters / 1000;
    colBounds = find(xCentersKm >= xLimits(1) & xCentersKm <= xLimits(2));
    rowBounds = find(yCentersKm >= yLimits(1) & yCentersKm <= yLimits(2));
    localTB = smosEase.tb(rowBounds, colBounds);

    fig = figure('Visible', 'off', 'Color', 'w', ...
        'Position', [80 80 1440 690]);
    tiles = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', ...
        'Padding', 'compact');

    ax1 = nexttile(tiles);
    h1 = imagesc(ax1, xCentersKm(colBounds), ...
        yCentersKm(rowBounds), localTB, [205 255]);
    h1.AlphaData = isfinite(localTB);
    set(ax1, 'YDir', 'normal', 'Color', [0.94 0.95 0.94]);
    hold(ax1, 'on');
    scatter(ax1, oibX, oibY, 5, [0.82 0.18 0.13], 'filled', ...
        'MarkerFaceAlpha', 0.55, 'MarkerEdgeAlpha', 0.55);
    axis(ax1, 'equal'); xlim(ax1, xLimits); ylim(ax1, yLimits);
    title(ax1, '① EASE 亮温与同日 OIB 航迹');
    xlabel(ax1, 'EASE x (km)'); ylabel(ax1, 'EASE y (km)');

    ax2 = nexttile(tiles);
    h2 = imagesc(ax2, xCentersKm(colBounds), ...
        yCentersKm(rowBounds), localTB, [205 255]);
    h2.AlphaData = isfinite(localTB);
    set(ax2, 'YDir', 'normal', 'Color', [0.94 0.95 0.94]);
    hold(ax2, 'on');
    scatter(ax2, oibX, oibY, 4, [0.72 0.72 0.72], 'filled', ...
        'MarkerFaceAlpha', 0.32, 'MarkerEdgeAlpha', 0.32);
    scatter(ax2, pairX, pairY, 28, 's', 'MarkerFaceColor', green, ...
        'MarkerEdgeColor', 'w', 'LineWidth', 0.45, ...
        'MarkerFaceAlpha', 0.78);
    axis(ax2, 'equal'); xlim(ax2, xLimits); ylim(ax2, yLimits);
    title(ax2, sprintf('② 3 × 3 汇总后的 %d 个配对中心格', height(pairs)));
    xlabel(ax2, 'EASE x (km)'); ylabel(ax2, 'EASE y (km)');

    colormap(fig, parula(256));
    cb = colorbar(ax2, 'eastoutside');
    cb.Label.String = 'SMOS TB (K)';
    set([ax1 ax2], 'FontName', fontName, 'FontSize', 11.5, ...
        'Layer', 'top', 'Box', 'on');
    title(tiles, '2012-03-14 EASE 网格上的同日配准结果', ...
        'FontName', fontName, ...
        'FontSize', 17, 'FontWeight', 'bold', 'Color', green);
    exportgraphics(fig, outputPath, 'Resolution', 190);
    close(fig);
end

function plot_daily_summary(summary, outputPath, fontName, green, grey, orange)
    excluded = summary.tbOverlapCells - summary.matchedCells;
    positions = (1:height(summary)).';
    displayDays = summary.day;
    displayDays.Format = 'yyyy-MM-dd';
    labels = string(displayDays);
    fig = figure('Visible', 'off', 'Color', 'w', ...
        'Position', [80 80 1450 620]);
    ax = axes(fig);
    bars = bar(ax, positions, [summary.matchedCells excluded], ...
        'stacked', 'BarWidth', 0.82);
    bars(1).FaceColor = green;
    bars(2).FaceColor = grey;
    hold(ax, 'on');
    noMatch = summary.matchedCells == 0;
    scatter(ax, positions(noMatch), 12 * ones(nnz(noMatch), 1), ...
        80, orange, 'x', 'LineWidth', 2.3);
    grid(ax, 'on');
    ylabel(ax, 'SMOS center cells');
    xlabel(ax, 'OIB flight day');
    title(ax, '23 个飞行日中 22 天形成有效配对，2012-03-28 无有效 SMOS 重叠');
    legend(ax, {'matched single-ice-type cells', ...
        'excluded mixed-type cells', 'no valid overlap'}, ...
        'Location', 'northwest');
    set(ax, 'FontName', fontName, 'FontSize', 11);
    xticks(ax, positions);
    xticklabels(ax, labels);
    xtickangle(ax, 45);
    xlim(ax, [0.2 height(summary) + 0.8]);
    exportgraphics(fig, outputPath, 'Resolution', 190);
    close(fig);
end

function plot_pair_overview(pairs, outputPath, fontName, green, lightGreen, orange, blue)
    years = year(pairs.day);
    fig = figure('Visible', 'off', 'Color', 'w', ...
        'Position', [80 80 1320 790]);
    tiles = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', ...
        'Padding', 'compact');

    ax1 = nexttile(tiles);
    histogram(ax1, pairs.TB(years == 2012), 20, ...
        'FaceColor', green, 'FaceAlpha', 0.72, 'EdgeColor', 'none');
    hold(ax1, 'on');
    histogram(ax1, pairs.TB(years == 2015), 20, ...
        'FaceColor', orange, 'FaceAlpha', 0.60, 'EdgeColor', 'none');
    xlabel(ax1, 'SMOS TB (K)'); ylabel(ax1, 'paired cells');
    title(ax1, '观测亮温分布');
    legend(ax1, {'2012', '2015'}, 'Location', 'northwest'); grid(ax1, 'on');

    ax2 = nexttile(tiles);
    histogram(ax2, log10(pairs.nOIB), 22, ...
        'FaceColor', blue, 'EdgeColor', 'none');
    xline(ax2, log10(median(pairs.nOIB)), '--', ...
        sprintf('median = %.0f', median(pairs.nOIB)), ...
        'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom');
    xlabel(ax2, 'log10(OIB samples per 3 × 3 window)');
    ylabel(ax2, 'paired cells'); title(ax2, '每个配对格的 OIB 样本量');
    grid(ax2, 'on');

    ax3 = nexttile(tiles);
    counts = [nnz(years == 2012 & ~pairs.isMYI), ...
        nnz(years == 2012 & pairs.isMYI); ...
        nnz(years == 2015 & ~pairs.isMYI), ...
        nnz(years == 2015 & pairs.isMYI)];
    bars = bar(ax3, categorical({'2012','2015'}), counts, 'stacked');
    bars(1).FaceColor = lightGreen; bars(2).FaceColor = green;
    ylabel(ax3, 'paired cells'); title(ax3, '单一冰型筛选后的构成');
    legend(ax3, {'FYI', 'MYI'}, 'Location', 'northwest'); grid(ax3, 'on');

    ax4 = nexttile(tiles);
    sample = 1:max(1, floor(height(pairs)/2600)):height(pairs);
    fyi = sample(~pairs.isMYI(sample));
    myi = sample(pairs.isMYI(sample));
    scatter(ax4, pairs.hs(fyi), pairs.hi(fyi), 10, lightGreen, ...
        'filled', 'MarkerFaceAlpha', 0.45); hold(ax4, 'on');
    scatter(ax4, pairs.hs(myi), pairs.hi(myi), 10, green, ...
        'filled', 'MarkerFaceAlpha', 0.45);
    xlabel(ax4, 'OIB mean snow depth h_s (m)');
    ylabel(ax4, 'OIB mean ice thickness h_i (m)');
    title(ax4, '配对格中的 OIB 验证量');
    legend(ax4, {'FYI', 'MYI'}, 'Location', 'northwest'); grid(ax4, 'on');

    set([ax1 ax2 ax3 ax4], 'FontName', fontName, 'FontSize', 10.5);
    title(tiles, sprintf('%d 条同日 OIB–SMOS 配对样本', height(pairs)), ...
        'FontName', fontName, 'FontSize', 17, 'FontWeight', 'bold');
    exportgraphics(fig, outputPath, 'Resolution', 190);
    close(fig);
end
