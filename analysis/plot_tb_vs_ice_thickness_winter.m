function result = plot_tb_vs_ice_thickness_winter( ...
        outputDirectory, limitBrinePermittivity)
%PLOT_TB_VS_ICE_THICKNESS_WINTER 复现论文图3.3的场景设计。
%   使用当前 tb_forward 的20层海冰、1层干雪配置，绘制FYI/MYI在
%   三种表面温度和三种雪厚下的平均L波段亮温随冰厚变化。
%
%   论文图3.3的场景：
%     FYI: h_s = 1, 10, 20 cm; h_i = 0--2 m
%     MYI: h_s = 1, 10, 25 cm; h_i = 0--4 m
%     T_surf = -15, -25, -35 degC
%
%   论文没有给出图3.3逐曲线的大气亮温修正，因此本图明确采用
%   atmosphereTBK=0 K。入射角为0:1:40度，H/V极化作算术平均。
%   LIMITBRINEPERMITTIVITY 默认 false，与冻结主线一致；只有显式传入
%   true 时才绘制把 Vb 截到 70 permille 的敏感性曲线。辐射传输始终
%   保持 Burke 模式。

    if nargin < 1 || strlength(string(outputDirectory)) == 0
        outputDirectory = fullfile(fileparts(mfilename('fullpath')), ...
            'tb_thickness_sensitivity');
    end
    if nargin < 2 || isempty(limitBrinePermittivity)
        limitBrinePermittivity = false;
    end
    if ~isfolder(outputDirectory)
        mkdir(outputDirectory);
    end

    projectDirectory = fileparts(fileparts(mfilename('fullpath')));
    addpath(projectDirectory);

    surfaceTemperaturesC = [-15, -25, -35];
    lineStyles = {':', '--', '-'};
    colors = [0.0000, 0.4470, 0.6980; ... % blue, Okabe-Ito
              0.0000, 0.6196, 0.4510; ... % bluish green
              0.8353, 0.3686, 0.0000];    % vermillion
    iceTypes = ["FYI", "MYI"];
    isMYI = [false, true];
    maxIceThicknessM = [2, 4];
    snowThicknessM = [0.01, 0.10, 0.20; 0.01, 0.10, 0.25];

    baseModel = struct( ...
        'surfaceTempC', -25, ...
        'waterTempC', -1.8, ...
        'waterSalinity', 33, ...
        'nIceLayers', 20, ...
        'nSnowLayers', 1, ...
        'incidenceAnglesDeg', 0:1:40, ...
        'includeCosmic', true, ...
        'cosmicTBK', 2.7, ...
        'atmosphereTBK', 0, ...
        'tbOffsetK', 0, ...
        'radiativeTransferMode', "burke", ...
        'limitBrinePermittivity', limitBrinePermittivity, ...
        'allowColdBrineExtrapolation', true);

    figureHandle = figure('Color', 'w', 'Units', 'centimeters', ...
        'Position', [2, 2, 19.0, 8.8]);
    layout = tiledlayout(figureHandle, 1, 2, 'TileSpacing', 'compact', ...
        'Padding', 'compact');

    nPoint = 201;
    nRows = 2 * 3 * 3 * nPoint;
    dataIceType = strings(nRows, 1);
    dataIsMYI = false(nRows, 1);
    dataIceThicknessM = zeros(nRows, 1);
    dataSnowThicknessM = zeros(nRows, 1);
    dataSurfaceTemperatureC = zeros(nRows, 1);
    dataBrightnessTemperatureK = zeros(nRows, 1);
    row = 1;

    axesHandles = gobjects(1, 2);
    for iType = 1:2
        axesHandle = nexttile(layout);
        axesHandles(iType) = axesHandle;
        hold(axesHandle, 'on');

        iceThicknessM = linspace(0.01, maxIceThicknessM(iType), nPoint);
        for iSnow = 1:3
            hs = snowThicknessM(iType, iSnow);
            for iTemp = 1:3
                surfaceTempC = surfaceTemperaturesC(iTemp);
                model = baseModel;
                model.surfaceTempC = surfaceTempC;
                tb = tb_forward(iceThicknessM, hs * ones(1, nPoint), ...
                    isMYI(iType), model);

                plot(axesHandle, iceThicknessM, tb, ...
                    'Color', colors(iSnow,:), ...
                    'LineStyle', lineStyles{iTemp}, ...
                    'LineWidth', 1.35);

                rows = row:(row + nPoint - 1);
                dataIceType(rows) = iceTypes(iType);
                dataIsMYI(rows) = isMYI(iType);
                dataIceThicknessM(rows) = iceThicknessM;
                dataSnowThicknessM(rows) = hs;
                dataSurfaceTemperatureC(rows) = surfaceTempC;
                dataBrightnessTemperatureK(rows) = tb;
                row = row + nPoint;
            end
        end

        xlim(axesHandle, [0, maxIceThicknessM(iType)]);
        ylim(axesHandle, [140, 260]);
        if isMYI(iType)
            xticks(axesHandle, 0:1:4);
        else
            xticks(axesHandle, 0:0.5:2);
        end
        yticks(axesHandle, 140:20:260);
        grid(axesHandle, 'on');
        axesHandle.GridAlpha = 0.15;
        axesHandle.MinorGridAlpha = 0.08;
        axesHandle.LineWidth = 0.8;
        axesHandle.FontName = 'Times New Roman';
        axesHandle.FontSize = 8.5;
        axesHandle.TickDir = 'out';
        box(axesHandle, 'on');
        xlabel(axesHandle, '海冰厚度  h_i  (m)', ...
            'FontName', 'Microsoft YaHei', 'FontSize', 9);
        title(axesHandle, sprintf('(%c) %s', ...
            char('a' + iType - 1), iceTypes(iType)), ...
            'FontName', 'Times New Roman', 'FontSize', 10, ...
            'FontWeight', 'normal');

        snowHandles = gobjects(1, 3);
        temperatureHandles = gobjects(1, 3);
        for index = 1:3
            snowHandles(index) = plot(axesHandle, NaN, NaN, '-', ...
                'Color', colors(index,:), 'LineWidth', 1.6);
            temperatureHandles(index) = plot(axesHandle, NaN, NaN, ...
                lineStyles{index}, 'Color', [0.15, 0.15, 0.15], ...
                'LineWidth', 1.4);
        end
        snowCm = round(100 * snowThicknessM(iType,:));
        legendLabels = [ ...
            compose('h_s = %d cm', snowCm), ...
            compose('T_{surf} = %d °C', surfaceTemperaturesC)];
        legend(axesHandle, [snowHandles, temperatureHandles], ...
            legendLabels, 'Location', 'southeast', 'NumColumns', 1, ...
            'Box', 'off', 'FontName', 'Times New Roman', ...
            'FontSize', 7.6, 'Interpreter', 'tex');
    end

    ylabel(axesHandles(1), '亮温  T_B  (K)', ...
        'FontName', 'Microsoft YaHei', 'FontSize', 9);
    if limitBrinePermittivity
        limitLabel = 'V_b上限70‰';
    else
        limitLabel = 'V_b不截断';
    end
    title(layout, sprintf( ...
        '典型北极冬季条件下亮温随海冰厚度的变化（Burke，%s）', ...
        limitLabel), ...
        'FontName', 'Microsoft YaHei', 'FontSize', 11, ...
        'FontWeight', 'normal');

    outputTable = table(dataIceType, dataIsMYI, dataIceThicknessM, ...
        dataSnowThicknessM, dataSurfaceTemperatureC, ...
        dataBrightnessTemperatureK, ...
        'VariableNames', {'iceType','isMYI','iceThicknessM', ...
        'snowThicknessM','surfaceTemperatureC','brightnessTemperatureK'});

    if limitBrinePermittivity
        baseName = 'tb_vs_ice_thickness_winter_current20layers';
    else
        baseName = 'tb_vs_ice_thickness_winter_burke_uncapped';
    end
    csvFile = fullfile(outputDirectory, [baseName, '.csv']);
    pngFile = fullfile(outputDirectory, [baseName, '.png']);
    figFile = fullfile(outputDirectory, [baseName, '.fig']);
    writetable(outputTable, csvFile);
    exportgraphics(figureHandle, pngFile, 'Resolution', 400);
    savefig(figureHandle, figFile);

    result = struct('figure', figureHandle, 'data', outputTable, ...
        'pngFile', string(pngFile), 'figFile', string(figFile), ...
        'csvFile', string(csvFile), 'model', baseModel);
    fprintf('冬季TB—冰厚曲线已保存：%s\n', pngFile);
end
