function output = plot_table1_tb_difference()
%PLOT_TABLE1_TB_DIFFERENCE 八个理想场景的冰厚/雪深—亮温差图。
%  只读取 processed/tb_forward_table1_validation.csv，绝不重算前向模型。
%  亮温差定义为 TBmodel - TBpaper；连线只辅助辨认同冰型离散场景。

    analysisRoot = fileparts(mfilename('fullpath'));
    matlabRoot = fileparts(analysisRoot);
    source = fullfile(matlabRoot, 'processed', ...
        'tb_forward_table1_validation.csv');
    data = readtable(source);
    required = {'scenario','isMYI','hi','hs','tbPaper','tbModel','differenceK'};
    assert(height(data) == 8 && all(ismember(required, data.Properties.VariableNames)), ...
        '需要论文表 1 的八个场景及完整亮温字段。');
    assert(numel(unique(string(data.scenario))) == 8 && ...
        all(isfinite(data.hi)) && all(isfinite(data.hs)) && ...
        all(isfinite(data.differenceK)) && ...
        all(abs(data.tbModel - data.tbPaper - data.differenceK) < 1e-9), ...
        '场景缺失或 differenceK 与 TBmodel-TBpaper 不一致。');

    colors = [0, 0.447, 0.698; 0.835, 0.369, 0];
    markers = {'o','s'};
    lineStyles = {'--',':'};
    labels = {'FYI（一年冰）','MYI（多年冰）'};
    xFields = {'hi','hs'};
    xLabels = {'冰厚 h_i (m)','雪深 h_s (m)'};
    fig = figure('Visible','off','Color','white', ...
        'Units','pixels','Position',[100,100,1200,600]);
    layout = tiledlayout(fig, 1, 2, 'TileSpacing','compact', ...
        'Padding','loose');
    for panel = 1:2
        ax = nexttile(layout, panel);
        hold(ax, 'on');
        lineHandles = gobjects(2, 1);
        for iceType = 0:1
            subset = data(data.isMYI == iceType, :);
            subset = sortrows(subset, xFields{panel});
            x = subset.(xFields{panel});
            y = subset.differenceK;
            lineHandles(iceType+1) = plot(ax, x, y, ...
                'Color',colors(iceType+1,:), ...
                'LineStyle',lineStyles{iceType+1}, 'LineWidth',1.3, ...
                'Marker',markers{iceType+1}, 'MarkerSize',7, ...
                'MarkerFaceColor',colors(iceType+1,:), ...
                'DisplayName',labels{iceType+1});
            xOffset = 0.010;
            if panel == 1, xOffset = 0.075; end
            for j = 1:height(subset)
                text(ax, x(j) + xOffset, y(j) + 0.08, ...
                    string(subset.scenario(j)), ...
                    'Color',colors(iceType+1,:), ...
                    'FontSize',10,'FontWeight','bold');
            end
        end
        yline(ax, 0, '-', 'Color',[0.55 0.55 0.55], ...
            'LineWidth',0.8, 'HandleVisibility','off');
        xlabel(ax, xLabels{panel});
        xlim(ax, [0, max(data.(xFields{panel})) * 1.12]);
        ylim(ax, [-4, 0.35]);
        grid(ax, 'on');
        ax.GridAlpha = 0.12;
        ax.FontName = 'Microsoft YaHei';
        ax.FontSize = 10;
        box(ax, 'off');
        if panel == 1
            ylabel(ax, '亮温差 TB_{model} - TB_{paper} (K)');
            legend(ax, lineHandles, labels, ...
                'Location','southwest','Box','off');
        end
    end
    sgtitle(layout, '8 个理想场景：冰厚、雪深与亮温差（离散场景连线）', ...
        'FontName','Microsoft YaHei','FontSize',15,'FontWeight','bold');

    stem = fullfile(analysisRoot, 'tb_difference_vs_hi_hs_8_scenarios');
    pngFile = [stem '.png'];
    svgFile = [stem '.svg'];
    %grayFile = [stem '_grayscale.png'];
    % Windows 的 MATLAB 图形后端有时不能直接覆盖已存在的 PNG。
    if isfile(pngFile), delete(pngFile); end
    if isfile(svgFile), delete(svgFile); end
    %if isfile(grayFile), delete(grayFile); end
    exportgraphics(fig, pngFile, 'Resolution',300);
    print(fig, svgFile, '-dsvg');
    % gray = uint8(0.299 * double(image(:,:,1)) + ...
    %     0.587 * double(image(:,:,2)) + 0.114 * double(image(:,:,3)));
    % imwrite(gray, grayFile);
    close(fig);
    output = struct('source',source,'png',pngFile,'svg',svgFile);
    %output = struct('source',source,'png',pngFile,'svg',svgFile, ...
    %    'grayscale',grayFile);
    fprintf('八场景图已保存：%s\n', pngFile);
end
