function plot_solution_space_2017(spaceFYI, spaceMYI, outputDir)
%PLOT_SOLUTION_SPACE_2017 绘制论文图 3/4 对应的两种输入解空间。
%   每张图按 FYI/MYI 与冰厚/雪深排成 2×2；使用前向模型计算的亮温。

    if ~isfolder(outputDir), mkdir(outputDir); end
    iceTypes = {spaceFYI, spaceMYI};
    branches = {'FBice', 'FBsnow'};
    fileNames = {'solution_space_radar_2017.png', ...
        'solution_space_laser_2017.png'};
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
                scatter(ax, space.TB, space.(branches{b}), 6, values, 'filled', ...
                    'MarkerEdgeColor', 'none');
                colormap(ax, parula);
                clim(ax, limits);
                cb = colorbar(ax);
                cb.Label.String = colorLabel;
                xlabel(ax, '模拟亮温 (K)');
                ylabel(ax, [branches{b} ' (m)']);
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
