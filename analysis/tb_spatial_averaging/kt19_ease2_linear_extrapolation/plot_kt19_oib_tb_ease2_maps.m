function plot_kt19_oib_tb_ease2_maps(pairs, selectedDays, pngPath)
%PLOT_KT19_OIB_TB_EASE2_MAPS 绘制代表日 OIB 模拟 TB 的 EASE2 分布。

    selectedDays = dateshift(datetime(selectedDays(:)),'start','day');
    selected = ismember(pairs.day,selectedDays) & ...
        isfinite(pairs.TBOIBKT19Linear);
    data = pairs(selected,:);
    if isempty(data)
        error('plot_kt19_oib_tb_ease2_maps:NoData', ...
            '所选日期没有有效的 OIB 模拟亮温。');
    end
    gridDefinition = grid_resampling.ease_grid("EASE2");
    values = data.TBOIBKT19Linear;
    colorLimits = [floor(min(values)), ceil(max(values))];
    if diff(colorLimits) < 1, colorLimits = colorLimits + [-0.5 0.5]; end
    days = unique(data.day).';
    figureHandle = figure('Visible','off','Color','w', ...
        'Position',[60 60 max(680,650*numel(days)) 570]);
    layout = tiledlayout(figureHandle,1,numel(days), ...
        'TileSpacing','compact','Padding','compact');
    title(layout,'OIB KT19 驱动的模拟亮温空间分布（EASE2-N）', ...
        'FontName','Microsoft YaHei','FontWeight','bold');
    for dayIndex = 1:numel(days)
        ax = nexttile(layout);
        current = data(data.day==days(dayIndex),:);
        xKm = gridDefinition.xCenters(current.col).' ./ 1000;
        yKm = gridDefinition.yCenters(current.row).' ./ 1000;
        scatter(ax,xKm,yKm,105,current.TBOIBKT19Linear,'s','filled', ...
            'MarkerEdgeColor','none');
        axis(ax,'equal'); box(ax,'on'); grid(ax,'on');
        colormap(ax,parula(256)); clim(ax,colorLimits);
        cb = colorbar(ax); cb.Label.String = '模拟亮温 (K)';
        xlabel(ax,'EASE2 x (km)'); ylabel(ax,'EASE2 y (km)');
        title(ax,sprintf('%s，n=%d\n平均表温 %.1f°C', ...
            char(string(days(dayIndex),'yyyy-MM-dd')),height(current), ...
            mean(current.surfaceTempC,'omitnan')), ...
            'FontName','Microsoft YaHei');
        ax.FontName='Microsoft YaHei'; ax.GridAlpha=0.18;
    end
    exportgraphics(figureHandle,pngPath,'Resolution',300);
    close(figureHandle);
end
