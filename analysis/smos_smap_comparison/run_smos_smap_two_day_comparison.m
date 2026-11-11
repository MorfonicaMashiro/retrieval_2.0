function [pairs, metrics, summary] = run_smos_smap_two_day_comparison()
%RUN_SMOS_SMAP_TWO_DAY_COMPARISON 比较 2015-04-01/03 的 SMOS 与 SMAP。
%   只使用当前实验统一口径的步长 3 OIB 覆盖中心；不运行前向或反演。
%   差值定义为 SMAP-SMOS。

    branchRoot = fileparts(mfilename('fullpath'));
    matlabRoot = fileparts(fileparts(branchRoot));
    workspaceRoot = fileparts(matlabRoot);
    dataRoot = fullfile(workspaceRoot, 'data');
    resultRoot = fullfile(branchRoot, 'results');
    figureRoot = fullfile(branchRoot, 'figures');
    addpath(matlabRoot, fullfile(matlabRoot, 'readdata'), ...
        fullfile(matlabRoot, 'plot'), branchRoot);
    if ~isfolder(resultRoot), mkdir(resultRoot); end
    if ~isfolder(figureRoot), mkdir(figureRoot); end

    targetDays = [datetime(2015,4,1); datetime(2015,4,3)];
    outputStem = "smos_smap_20150401_20150403";
    operatorOptions = struct('centerStride', 3, 'centerPhase', [2 2]);
    pairs = table();
    summary = table();

    for dayIndex = 1:numel(targetDays)
        day = targetDays(dayIndex);
        token = char(string(day, 'yyyyMMdd'));
        oibFile = fullfile(dataRoot, 'OIB', 'NSIDC-0708_QuickLook', ...
            '2015_GR_NASA', ['OIB_' token '_IDCSI2.txt']);
        smosFile = fullfile(dataRoot, 'smos', '2015', ...
            ['SMOS_TB_north_' token '.nc']);
        smapDirectory = fullfile(dataRoot, 'smap', 'SPL1CTB_006', ...
            '2015', char(string(day,'MM')), char(string(day,'dd')));
        require_file(oibFile, 'OIB');
        require_file(smosFile, 'SMOS');
        smapListing = dir(fullfile(smapDirectory, ...
            ['SMAP_L1C_TB_*_' token 'T*_R19240_001.h5']));
        if isempty(smapListing)
            error('run_smos_smap_two_day_comparison:MissingSMAP', ...
                ['%s 没有找到 SMAP HDF5。请按 smap_download_manifest.csv ' ...
                '下载到：%s'], token, smapDirectory);
        end
        smapFiles = string(fullfile({smapListing.folder}, ...
            {smapListing.name})).';

        oib = read_oib(oibFile);
        smosEase = grid_resampling.regrid_smos_to_ease( ...
            read_smos(smosFile), "EASE");
        [dailyPairs, collocateStats] = collocate_day( ...
            oib, smosEase, day, [], [], operatorOptions);
        smapDaily = build_smap_daily_composite(smapFiles);
        dailyPairs = attach_smap(dailyPairs, smosEase.grid, smapDaily);
        pairs = [pairs; dailyPairs]; %#ok<AGROW>

        valid = isfinite(dailyPairs.TB) & isfinite(dailyPairs.TBSMAPMeanHV);
        summary = [summary; table(day, numel(smapFiles), ...
            collocateStats.rawPoints, collocateStats.qualityPoints, ...
            height(dailyPairs), nnz(valid), nnz(~valid), ...
            sum(dailyPairs.nSMAPLooks(valid)), ...
            'VariableNames', {'day','nSMAPFiles','rawOIBPoints', ...
            'qualityOIBPoints','stride3SMOSPairs','commonPairs', ...
            'missingSMAPPairs','totalSMAPLooksAtCommonPairs'})]; %#ok<AGROW>
        fprintf('%s：SMAP 半轨 %d，步长3 SMOS 配对 %d，共同有效 %d。\n', ...
            token, numel(smapFiles), height(dailyPairs), nnz(valid));
    end

    pairs = sortrows(pairs, {'day','row','col'});
    metrics = build_metrics(pairs);
    matPath = fullfile(resultRoot, outputStem + ".mat");
    save(matPath, 'pairs', 'metrics', 'summary', 'operatorOptions', '-v7.3');
    writetable(pairs, fullfile(resultRoot, outputStem + "_pairs.csv"));
    writetable(metrics, fullfile(resultRoot, outputStem + "_metrics.csv"));
    writetable(summary, fullfile(resultRoot, outputStem + "_summary.csv"));
    plot_comparison(pairs, fullfile(figureRoot, ...
        outputStem + "_comparison.png"));
    fprintf('SMOS-SMAP 两日对比完成：%d 条共同有效配对。\n%s\n', ...
        nnz(isfinite(pairs.TB) & isfinite(pairs.TBSMAPMeanHV)), matPath);
end

function pairs = attach_smap(pairs, smosGrid, smap)
    xOld = smosGrid.xCenters(double(pairs.col));
    yOld = smosGrid.yCenters(double(pairs.row));
    xOld = xOld(:);
    yOld = yOld(:);
    [latitude, longitude] = projinv(smosGrid.projection, xOld, yOld);
    [x, y] = projfwd(smap.projection, latitude, longitude);
    smapRow = floor((smap.edge - y) ./ smap.cellSize) + 1;
    smapCol = floor((x + smap.edge) ./ smap.cellSize) + 1;
    inside = smapRow >= 1 & smapRow <= smap.nRows & ...
        smapCol >= 1 & smapCol <= smap.nCols & ...
        isfinite(smapRow) & isfinite(smapCol);
    linear = nan(height(pairs),1);
    linear(inside) = sub2ind([smap.nRows smap.nCols], ...
        smapRow(inside), smapCol(inside));
    tbH = nan(height(pairs),1);
    tbV = nan(height(pairs),1);
    tbMeanHV = nan(height(pairs),1);
    nLooks = zeros(height(pairs),1);
    tbH(inside) = smap.tbH(linear(inside));
    tbV(inside) = smap.tbV(linear(inside));
    tbMeanHV(inside) = smap.tbMeanHV(linear(inside));
    nLooks(inside) = smap.observationCount(linear(inside));
    smapCenterX = nan(height(pairs),1);
    smapCenterY = nan(height(pairs),1);
    smapCenterX(inside) = smap.xCenters(smapCol(inside));
    smapCenterY(inside) = smap.yCenters(smapRow(inside));
    centerDistanceKm = hypot(x - smapCenterX, y - smapCenterY) / 1000;
    maximumDistanceKm = sqrt(2) * smap.cellSize / 2 / 1000 + 1e-6;
    tooFar = centerDistanceKm > maximumDistanceKm;
    tbH(tooFar) = NaN; tbV(tooFar) = NaN; tbMeanHV(tooFar) = NaN;
    nLooks(tooFar) = 0;
    difference = tbMeanHV - pairs.TB;
    meanTB = (tbMeanHV + pairs.TB) / 2;
    pairs = addvars(pairs, latitude, longitude, smapRow, smapCol, ...
        centerDistanceKm, nLooks, tbH, tbV, tbMeanHV, difference, meanTB, ...
        'NewVariableNames', {'latitude','longitude','smapRow','smapCol', ...
        'smapCenterDistanceKm','nSMAPLooks','TBSMAPH','TBSMAPV', ...
        'TBSMAPMeanHV','TBSMAPMinusSMOS','TBPairMean'});
end

function metrics = build_metrics(pairs)
    days = [NaT; unique(pairs.day)];
    channels = ["H"; "V"; "(H+V)/2"];
    values = {pairs.TBSMAPH, pairs.TBSMAPV, pairs.TBSMAPMeanHV};
    metrics = table();
    for dayIndex = 1:numel(days)
        if ismissing(days(dayIndex))
            dayMask = true(height(pairs),1);
            label = "两日合并";
        else
            dayMask = pairs.day == days(dayIndex);
            label = string(days(dayIndex), 'yyyy-MM-dd');
        end
        for channelIndex = 1:numel(channels)
            x = pairs.TB(dayMask);
            y = values{channelIndex}(dayMask);
            fit = fit_retrieval(x, y);
            valid = isfinite(x) & isfinite(y);
            differenceStd = std(y(valid) - x(valid), 0, 'omitnan');
            current = table(label, channels(channelIndex), fit.n, ...
                fit.correlation, fit.rSquared, fit.bias, fit.rmse, ...
                differenceStd, fit.slope, fit.intercept, ...
                'VariableNames', {'period','SMAPChannel','n', ...
                'correlation','rSquared','biasSMAPMinusSMOS','rmse', ...
                'differenceStd','slope','intercept'});
            metrics = [metrics; current]; %#ok<AGROW>
        end
    end
end

function plot_comparison(pairs, pngPath)
    valid = isfinite(pairs.TB) & isfinite(pairs.TBSMAPMeanHV);
    if ~any(valid)
        error('run_smos_smap_two_day_comparison:NoCommonData', ...
            '没有共同有效的 SMOS-SMAP 配对。');
    end
    colors = [0 114 189; 230 126 0] / 255;
    markers = {'o','^'};
    days = unique(pairs.day(valid));
    fig = figure('Visible','off','Color','w','Position',[80 80 1540 500]);
    layout = tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');
    title(layout,'2015 年两日 SMOS–SMAP 观测亮温对比（步长 3）', ...
        'FontName','Microsoft YaHei','FontWeight','bold');

    ax1 = nexttile(layout); hold(ax1,'on');
    ax2 = nexttile(layout); hold(ax2,'on');
    ax3 = nexttile(layout); hold(ax3,'on');
    for dayIndex = 1:numel(days)
        selected = valid & pairs.day == days(dayIndex);
        label = char(string(days(dayIndex),'yyyy-MM-dd'));
        scatter(ax1,pairs.TB(selected),pairs.TBSMAPMeanHV(selected),34, ...
            colors(dayIndex,:),markers{dayIndex},'filled', ...
            'MarkerFaceAlpha',0.72,'DisplayName',label);
        scatter(ax2,pairs.TBPairMean(selected), ...
            pairs.TBSMAPMinusSMOS(selected),34,colors(dayIndex,:), ...
            markers{dayIndex},'filled','MarkerFaceAlpha',0.72, ...
            'DisplayName',label);
        xJitter = dayIndex + linspace(-0.12,0.12,nnz(selected)).';
        scatter(ax3,xJitter,pairs.TBSMAPMinusSMOS(selected),34, ...
            colors(dayIndex,:),markers{dayIndex},'filled', ...
            'MarkerFaceAlpha',0.72,'DisplayName',label);
        medianDifference = median(pairs.TBSMAPMinusSMOS(selected),'omitnan');
        plot(ax3,[dayIndex-0.22 dayIndex+0.22], ...
            [medianDifference medianDifference],'-','Color',colors(dayIndex,:), ...
            'LineWidth',2.2,'HandleVisibility','off');
    end

    bounds = [min([pairs.TB(valid);pairs.TBSMAPMeanHV(valid)]), ...
        max([pairs.TB(valid);pairs.TBSMAPMeanHV(valid)])];
    plot(ax1,bounds,bounds,'k-','LineWidth',1.1,'DisplayName','1:1');
    xlim(ax1,bounds); ylim(ax1,bounds); axis(ax1,'square');
    fit = fit_retrieval(pairs.TB(valid),pairs.TBSMAPMeanHV(valid));
    text(ax1,0.04,0.96,sprintf(['n=%d  R^2=%.2f\n' ...
        'Bias=%+.2f K  RMSE=%.2f K'], ...
        fit.n,fit.rSquared,fit.bias,fit.rmse),'Units','normalized', ...
        'VerticalAlignment','top','BackgroundColor','w','Margin',4, ...
        'FontName','Microsoft YaHei');
    xlabel(ax1,'SMOS 日亮温 (K)'); ylabel(ax1,'SMAP 日均 (H+V)/2 (K)');
    title(ax1,'共同样本散点'); legend(ax1,'Location','southeast');

    difference = pairs.TBSMAPMinusSMOS(valid);
    bias = mean(difference,'omitnan');
    spread = std(difference,0,'omitnan');
    yline(ax2,0,'k-','LineWidth',1.0,'HandleVisibility','off');
    yline(ax2,bias,'-','Color',[0.00 0.45 0.70],'LineWidth',1.8, ...
        'DisplayName',sprintf('均值 %+.2f K',bias));
    yline(ax2,bias + 1.96*spread,'--','Color',[0.80 0.33 0.00], ...
        'LineWidth',1.2,'DisplayName','均值 ± 1.96 SD');
    yline(ax2,bias - 1.96*spread,'--','Color',[0.80 0.33 0.00], ...
        'LineWidth',1.2,'HandleVisibility','off');
    xlabel(ax2,'两传感器均值 (K)'); ylabel(ax2,'SMAP - SMOS (K)');
    title(ax2,'差值与亮温水平'); legend(ax2,'Location','best');

    yline(ax3,0,'k-','LineWidth',1.0,'HandleVisibility','off');
    xlim(ax3,[0.5 numel(days)+0.5]);
    xticks(ax3,1:numel(days));
    xticklabels(ax3,cellstr(string(days,'MM-dd')));
    ylabel(ax3,'SMAP - SMOS (K)'); title(ax3,'逐日差值分布');

    axesList = [ax1 ax2 ax3];
    for ax = axesList
        ax.FontName = 'Microsoft YaHei';
        ax.GridColor = [0.82 0.90 0.95];
        ax.GridAlpha = 0.65;
        grid(ax,'on'); box(ax,'on');
    end
    exportgraphics(fig,pngPath,'Resolution',300);
    close(fig);
end

function require_file(filePath, label)
    if ~isfile(filePath)
        error('run_smos_smap_two_day_comparison:MissingInput', ...
            '%s file not found: %s', label, filePath);
    end
end
