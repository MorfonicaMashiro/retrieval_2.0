function plot_2017_pairs(pairs, gridDefinition, outputPath)
%PLOT_2017_PAIRS 只绘制配对位置与每格 OIB 样本数，供空间配准检查。

    if isempty(pairs), error('plot_2017_pairs:NoData', '没有有效配对。'); end
    fig = figure('Visible', 'off', 'Color', 'w', ...
        'Position', [100 100 1050 570]);
    tiles = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact');
    years = [2012 2015];
    for k = 1:2
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
