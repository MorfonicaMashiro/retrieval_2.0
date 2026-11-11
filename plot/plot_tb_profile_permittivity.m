%% PLOT_TB_PROFILE_PERMITTIVITY
% 绘制前向亮温模型的两类中间量：
%   1) 固定几何与边界条件下的 20 层温度、盐度和海冰复介电常数廓线；
%   2) 温度、盐度到卤水体积的二维响应，以及 Vb 到介电常数的线性映射。
%
% 受控示例参数：hi=1.5 m、hs=0.15 m、Ts=-30 degC、Tw=-1.8 degC。
% 其中 MYI 几何对应 Data Synergy (2017) Table 1 场景 V；FYI 使用相同
% 几何用于控制变量比较，不代表额外观测场景。
%
% 公式来源：
%   温度廓线：周璐博士论文 Sect. 3.2.2, Eqs. (3-1)--(3-4)
%   FYI/MYI 盐度：同上 Eqs. (3-5)--(3-6)
%   卤水体积：同上 Eqs. (3-10)--(3-13)，Cox and Weeks (1983)
%   海冰介电常数：同上 Eq. (3-9)，Vant et al. (1978)

clear;
close all;

scriptDir = fileparts(mfilename('fullpath'));
projectDir = fileparts(scriptDir);
outputDir = fullfile(projectDir, 'figures', 'group_report');
if ~isfolder(outputDir)
    mkdir(outputDir);
end

set(groot, 'defaultAxesFontName', 'Microsoft YaHei');
set(groot, 'defaultTextFontName', 'Microsoft YaHei');
set(groot, 'defaultAxesFontSize', 11);
set(groot, 'defaultLineLineWidth', 1.8);

% 与前向模型默认值一致。
hi = 1.5;
hs = 0.15;
surfaceTempC = -30;
waterTempC = -1.8;
nIce = 20;
nSnow = 1;
snowConductivity = 0.31;
fyiMinSalinityPpt = 1.5;

fyi = build_profile(false, hi, hs, surfaceTempC, waterTempC, ...
    nIce, nSnow, snowConductivity, fyiMinSalinityPpt);
myi = build_profile(true, hi, hs, surfaceTempC, waterTempC, ...
    nIce, nSnow, snowConductivity, fyiMinSalinityPpt);

%% 图 1：模型实际使用的垂向输入与派生介电常数
fig1 = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [100, 100, 1600, 820]);
layout1 = tiledlayout(fig1, 2, 2, 'TileSpacing', 'compact', ...
    'Padding', 'compact');

colorFYI = [0.10, 0.45, 0.70];
colorMYI = [0.85, 0.37, 0.09];
colorSnow = [0.86, 0.93, 0.96];

ax = nexttile(layout1, 1);
hold(ax, 'on');
patch(ax, [-35, 0, 0, -35], [0, 0, hs, hs], colorSnow, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.65, 'HandleVisibility', 'off');
plot(ax, [surfaceTempC, fyi.snowTemperatureC, fyi.iceTemperatureC, waterTempC], ...
    [0, fyi.snowDepth, fyi.iceDepth, hs + hi], '-o', ...
    'Color', colorFYI, 'MarkerFaceColor', 'w', 'MarkerSize', 4, ...
    'DisplayName', 'FYI');
plot(ax, [surfaceTempC, myi.snowTemperatureC, myi.iceTemperatureC, waterTempC], ...
    [0, myi.snowDepth, myi.iceDepth, hs + hi], '--s', ...
    'Color', colorMYI, 'MarkerFaceColor', 'w', 'MarkerSize', 4, ...
    'DisplayName', 'MYI');
yline(ax, hs, ':', '雪/冰界面', 'Color', [0.35, 0.35, 0.35], ...
    'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
set(ax, 'YDir', 'reverse');
xlabel(ax, '温度 (°C)');
ylabel(ax, '距雪面的深度 (m)');
title(ax, '(a) 稳态温度廓线');
legend(ax, 'Location', 'southeast', 'Box', 'off');
grid(ax, 'on');

ax = nexttile(layout1, 2);
hold(ax, 'on');
plot(ax, fyi.salinityPpt, fyi.normalizedIceDepth, '-o', ...
    'Color', colorFYI, 'MarkerFaceColor', 'w', 'MarkerSize', 4, ...
    'DisplayName', 'FYI');
plot(ax, myi.salinityPpt, myi.normalizedIceDepth, '--s', ...
    'Color', colorMYI, 'MarkerFaceColor', 'w', 'MarkerSize', 4, ...
    'DisplayName', 'MYI');
set(ax, 'YDir', 'reverse');
xlabel(ax, '海冰盐度 (g kg^{-1})');
ylabel(ax, '归一化冰深 z/h_i');
title(ax, '(b) 海冰盐度廓线');
legend(ax, 'Location', 'southeast', 'Box', 'off');
grid(ax, 'on');

ax = nexttile(layout1, 3);
hold(ax, 'on');
plot(ax, real(fyi.epsIce), fyi.normalizedIceDepth, '-o', ...
    'Color', colorFYI, 'MarkerFaceColor', 'w', 'MarkerSize', 4, ...
    'DisplayName', 'FYI');
plot(ax, real(myi.epsIce), myi.normalizedIceDepth, '--s', ...
    'Color', colorMYI, 'MarkerFaceColor', 'w', 'MarkerSize', 4, ...
    'DisplayName', 'MYI');
set(ax, 'YDir', 'reverse');
xlabel(ax, '介电常数实部 ε''');
ylabel(ax, '归一化冰深 z/h_i');
title(ax, '(c) ε'' 垂向变化');
legend(ax, 'Location', 'southeast', 'Box', 'off');
grid(ax, 'on');

ax = nexttile(layout1, 4);
hold(ax, 'on');
plot(ax, imag(fyi.epsIce), fyi.normalizedIceDepth, '-o', ...
    'Color', colorFYI, 'MarkerFaceColor', 'w', 'MarkerSize', 4, ...
    'DisplayName', 'FYI');
plot(ax, imag(myi.epsIce), myi.normalizedIceDepth, '--s', ...
    'Color', colorMYI, 'MarkerFaceColor', 'w', 'MarkerSize', 4, ...
    'DisplayName', 'MYI');
set(ax, 'YDir', 'reverse');
xlabel(ax, '介电常数虚部 ε''''');
ylabel(ax, '归一化冰深 z/h_i');
title(ax, '(d) ε'''' 垂向变化');
legend(ax, 'Location', 'southeast', 'Box', 'off');
grid(ax, 'on');

sgtitle(layout1, sprintf(['20 层受控示例：h_i=%.1f m，h_s=%.2f m，' ...
    'T_s=%g °C，T_w=%g °C'], hi, hs, surfaceTempC, waterTempC), ...
    'FontWeight', 'bold', 'FontSize', 15);

exportgraphics(fig1, fullfile(outputDir, 'temperature_salinity_profiles.png'), ...
    'Resolution', 300, 'BackgroundColor', 'white');
close(fig1);

%% 图 2：复介电常数对温度和盐度的二维敏感性
temperatureAxis = linspace(-30, -0.25, 360);
salinityAxis = linspace(0.5, 8, 260);
[temperatureGrid, salinityGrid] = meshgrid(temperatureAxis, salinityAxis);
vbGrid = brine_volume_permille(temperatureGrid, salinityGrid);
valid = isfinite(vbGrid) & vbGrid >= 0 & vbGrid < 70;

vbShown = vbGrid;
vbShown(~valid) = NaN;
vbAxis = linspace(0, 70, 300);
epsRealLine = 3.10 + 0.0084 .* vbAxis;
epsImagFYI = 0.037 + 0.00445 .* vbAxis;
epsImagMYI = 0.003 + 0.00435 .* vbAxis;

fig2 = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [100, 100, 1600, 820]);
layout2 = tiledlayout(fig2, 2, 3, 'TileSpacing', 'compact', ...
    'Padding', 'compact');

plot_vb_heatmap(nexttile(layout2, 1, [2, 2]), temperatureAxis, ...
    salinityAxis, vbShown, fyi, myi);

ax = nexttile(layout2, 3);
plot(ax, vbAxis, epsRealLine, '-', 'Color', [0.10, 0.45, 0.70], ...
    'DisplayName', 'ε'' = 3.10 + 0.0084 V_b');
xlabel(ax, '卤水体积 V_b (‰)');
ylabel(ax, '介电常数实部 ε''');
title(ax, '(b) 实部：FYI 与 MYI 共用同一关系');
xlim(ax, [0, 70]);
grid(ax, 'on');
legend(ax, 'Location', 'northwest', 'Box', 'off');

ax = nexttile(layout2, 6);
hold(ax, 'on');
plot(ax, vbAxis, epsImagFYI, '-', 'Color', [0.10, 0.45, 0.70], ...
    'DisplayName', 'FYI: 0.037 + 0.00445 V_b');
plot(ax, vbAxis, epsImagMYI, '--', 'Color', [0.85, 0.37, 0.09], ...
    'DisplayName', 'MYI: 0.003 + 0.00435 V_b');
xlabel(ax, '卤水体积 V_b (‰)');
ylabel(ax, '介电常数虚部 ε''''');
title(ax, '(c) 虚部：冰型仅改变线性系数');
xlim(ax, [0, 70]);
grid(ax, 'on');
legend(ax, 'Location', 'northwest', 'Box', 'off');

sgtitle(layout2, ['温度和盐度先控制卤水体积，' ...
    '介电常数随后随 V_b 线性变化'], ...
    'FontWeight', 'bold', 'FontSize', 15);

exportgraphics(fig2, fullfile(outputDir, 'dielectric_sensitivity.png'), ...
    'Resolution', 300, 'BackgroundColor', 'white');
close(fig2);

fprintf('Saved profile and permittivity figures to %s\n', outputDir);

%% Local functions
function profile = build_profile(isMYI, hi, hs, surfaceTempC, waterTempC, ...
        nIce, nSnow, snowConductivity, fyiMinSalinityPpt)
    normalizedIceDepth = ((1:nIce) - 0.5) / nIce;
    if isMYI
        salinityPpt = 0.5 * 3.2 .* ...
            (1 - cos(pi .* normalizedIceDepth .^ ...
            (0.407 ./ (normalizedIceDepth + 0.573))));
    else
        bulkSalinity = 6.08 * exp(-5.81 * hi) + ...
            7.409 * exp(-0.5228 * hi);
        bulkSalinity = max(fyiMinSalinityPpt, bulkSalinity);
        salinityPpt = repmat(bulkSalinity, 1, nIce);
    end

    meanColumnTempC = 0.5 * (surfaceTempC + waterTempC);
    kIce = 2.034 + 0.13 * mean(salinityPpt) / meanColumnTempC;
    totalResistance = hs / snowConductivity + hi / kIce;
    heatFlux = (waterTempC - surfaceTempC) / totalResistance;
    interfaceTempC = surfaceTempC + heatFlux * hs / snowConductivity;

    snowDepth = (((1:nSnow) - 0.5) / nSnow) * hs;
    snowTemperatureC = surfaceTempC + ...
        heatFlux .* snowDepth / snowConductivity;
    iceDepth = hs + normalizedIceDepth * hi;
    iceTemperatureC = interfaceTempC + ...
        heatFlux .* (normalizedIceDepth * hi) / kIce;

    vbPermille = brine_volume_permille(iceTemperatureC, salinityPpt);
    if isMYI
        epsIce = 3.10 + 0.0084 .* vbPermille + ...
            1i .* (0.003 + 0.00435 .* vbPermille);
    else
        epsIce = 3.10 + 0.0084 .* vbPermille + ...
            1i .* (0.037 + 0.00445 .* vbPermille);
    end

    profile = struct('snowDepth', snowDepth, ...
        'snowTemperatureC', snowTemperatureC, ...
        'iceDepth', iceDepth, ...
        'normalizedIceDepth', normalizedIceDepth, ...
        'iceTemperatureC', iceTemperatureC, ...
        'salinityPpt', salinityPpt, ...
        'vbPermille', vbPermille, ...
        'epsIce', epsIce);
end

function vbPermille = brine_volume_permille(temperatureC, salinityPpt)
    f1 = nan(size(temperatureC));
    warm = temperatureC > -2;
    middle = temperatureC <= -2 & temperatureC > -22.9;
    cold = temperatureC <= -22.9;

    f1(warm) = polyval([0.21454, 0.58402, -18.407, -0.041221], ...
        temperatureC(warm));
    f1(middle) = polyval([-0.01074, -0.6397, -22.45, -4.732], ...
        temperatureC(middle));
    f1(cold) = polyval([0.7160, 55.27, 1309, 9899], ...
        temperatureC(cold));

    pureIceDensity = 0.917 - 1.403e-4 .* temperatureC;
    vbPermille = 1000 .* pureIceDensity .* salinityPpt ./ f1;
    % 旧的 F2 零气体体积闭合保留作追溯，不用于论文默认口径：
    % f2 = nan(size(temperatureC));
    % f2(warm) = polyval([1.3603e-4, 1.2291e-4, -0.016111, 0.090312], ...
    %     temperatureC(warm));
    % f2(middle) = polyval([-8.801e-6, -5.330e-4, -0.01763, 0.08903], ...
    %     temperatureC(middle));
    % f2(cold) = polyval([5.819e-4, 0.04518, 1.089, 8.547], ...
    %     temperatureC(cold));
    % vbPermille = 1000 .* pureIceDensity .* salinityPpt ./ ...
    %     (f1 - pureIceDensity .* salinityPpt .* f2);
end

function plot_vb_heatmap(ax, temperatureAxis, salinityAxis, values, fyi, myi)
    imageHandle = imagesc(ax, temperatureAxis, salinityAxis, values);
    imageHandle.AlphaData = ~isnan(values);
    set(ax, 'YDir', 'normal');
    set(ax, 'Color', 'w');
    colormap(ax, parula(256));
    colorBar = colorbar(ax);
    colorBar.Label.String = 'V_b (‰)';
    hold(ax, 'on');
    plot(ax, fyi.iceTemperatureC, fyi.salinityPpt, 'o', ...
        'Color', [0.10, 0.45, 0.70], 'MarkerFaceColor', 'w', ...
        'MarkerSize', 4.5, 'LineWidth', 1.2, 'DisplayName', 'FYI 20层状态');
    plot(ax, myi.iceTemperatureC, myi.salinityPpt, '^', ...
        'Color', [0.85, 0.37, 0.09], 'MarkerFaceColor', 'w', ...
        'MarkerSize', 4.5, 'LineWidth', 1.2, 'DisplayName', 'MYI 20层状态');
    xline(ax, -22.9, ':', 'Color', [0.25, 0.25, 0.25], ...
        'HandleVisibility', 'off');
    xline(ax, -2, ':', 'Color', [0.25, 0.25, 0.25], ...
        'HandleVisibility', 'off');
    xlabel(ax, '冰温 (°C)');
    ylabel(ax, '海冰盐度 (g kg^{-1})');
    title(ax, '(a) 非线性来源：V_b(T_i, S_i)，白色区域为 V_b ≥ 70‰');
    legend(ax, 'Location', 'northwest', 'Box', 'off');
end
