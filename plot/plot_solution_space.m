function plot_solution_space(spaceFYI, spaceMYI, outputDir)
%PLOT_SOLUTION_SPACE 绘制两种干舷输入对应的解空间。
%   每张图按 FYI/MYI 与冰厚/雪深排成 2×2；使用前向模型计算的亮温。

    if ~isfolder(outputDir), mkdir(outputDir); end
    iceTypes = {spaceFYI, spaceMYI};
    branches = {'FBice', 'FBsnow'};
    fileNames = {'solution_space_radar.png', ...
        'solution_space_laser.png'};
    for b = 1:2
        fig = figure('Visible', 'off', 'Color', 'w', ...
            'Position', [100 100 1100 780]);
        layout = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', ...
            'Padding', 'compact');
        for iceType = 1:2
            space = iceTypes{iceType};
            for variable = 1:2
                ax = nexttile(layout);
                if variable == 1
                    values = space.hi;
                    colorLabel = '冰厚 (m)';
                    limits = [0 10];
                else
                    values = space.hs;
                    colorLabel = '雪深 (m)';
                    limits = [0 1];
                end
                scatter(ax, space.(branches{b}), space.TB, 6, values, 'filled', ...
                    'MarkerEdgeColor', 'none');
                colormap(ax, parula);
                clim(ax, limits);
                cb = colorbar(ax);
                cb.Label.String = colorLabel;
                xlabel(ax, [branches{b} ' (m)']);
                ylabel(ax, '模拟亮温 (K)');
                if iceType == 1, typeName = 'FYI'; else, typeName = 'MYI'; end
                title(ax, [typeName '：' colorLabel]);
                grid(ax, 'on');
                ax.Box = 'on';
            end
        end
        exportgraphics(fig, fullfile(outputDir, fileNames{b}), 'Resolution', 160);
        close(fig);
    end
end
