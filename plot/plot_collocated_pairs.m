function plot_collocated_pairs(pairs, gridDefinition, outputPath)
%PLOT_COLLOCATED_PAIRS 绘制配对位置与每格 OIB 样本数。

    if isempty(pairs), error('plot_collocated_pairs:NoData', '没有有效配对。'); end
    years = unique(year(pairs.day)).';
    nColumns = min(2, numel(years));
    nRows = ceil(numel(years) / nColumns);
    fig = figure('Visible', 'off', 'Color', 'w', ...
        'Position', [100 100 1050 470 * nRows]);
    tiles = tiledlayout(fig, nRows, nColumns, 'TileSpacing', 'compact', ...
        'Padding', 'compact');
    for k = 1:numel(years)
        subset = year(pairs.day) == years(k);
        ax = nexttile(tiles);
        x = gridDefinition.xCenters(pairs.col(subset)) / 1000;
        y = gridDefinition.yCenters(pairs.row(subset)) / 1000;
        scatter(ax, x, y, 9, log10(pairs.nOIB(subset)), 'filled', ...
            'MarkerEdgeColor', 'none');
        axis(ax, 'equal');
        grid(ax, 'on');
        colormap(ax, parula);
        clim(ax, [0 4]);
        xlabel(ax, 'EASE x (km)');
        ylabel(ax, 'EASE y (km)');
        title(ax, sprintf('%d: %d flight days, %d cells', years(k), ...
            numel(unique(pairs.day(subset))), nnz(subset)));
    end
    cb = colorbar(ax);
    cb.Label.String = 'log10(OIB samples per 3x3 area)';
    outputPath = char(string(outputPath));
    if ~isfolder(fileparts(outputPath)), mkdir(fileparts(outputPath)); end
    exportgraphics(fig, outputPath, 'Resolution', 150);
    close(fig);
end
