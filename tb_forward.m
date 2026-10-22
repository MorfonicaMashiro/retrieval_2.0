function [tb, diagnostics] = tb_forward(hi, hs, isMYI, model)
%TB_FORWARD L 波段雪-冰-水前向亮温模型。
%   TB = TB_FORWARD(HI, HS, ISMYI, MODEL) 计算 1.4 GHz 亮温（K）。
%   HI、HS 单位为 m；ISMYI 为 0（FYI）或 1（MYI）。HI 和 HS 可以是
%   同形数组，ISMYI 当前为标量。默认结果是 0:1:40 度、H/V 两种极化的
%   算术平均，与 Zhou 等（2017）和 SMOS 日平均产品的约定一致。
%
%   默认物理配置复现 Zhou 等改进模型的公开描述：
%     * 海冰 20 层、干雪 1 层；
%     * 冰雪内稳态线性温度廓线；
%     * FYI 厚度相关体盐度，MYI Schwarzacher 型垂直盐度廓线；
%     * Vant 等海冰介电关系、Tiuri 等干雪介电关系；
%     * Klein-Swift 海水介电关系；
%     * Burke 等一阶非相干分层辐射传输递推。
%
%   MODEL 必需字段：surfaceTempC、waterTempC、waterSalinity。可选字段：
%     nIceLayers          默认 20；设为 1 可作单层敏感性基线
%     nSnowLayers         默认 1（博士论文明确采用 1 层雪）
%     frequencyHz         默认 1.4e9
%     incidenceAnglesDeg  默认 0:1:40
%     includeCosmic       默认 true
%     cosmicTBK           默认 2.7
%     atmosphereTBK       默认 0；原文未给出逐场景确定值
%     tbOffsetK           默认 0；不把验证偏差当作物理修正项
%     snowDensityKgM3     默认 320
%     snowConductivity    默认 0.31 W m^-1 K^-1
%     fyiMinSalinityPpt   默认 1.5
%     maxBrinePermille    默认 70（Vant 经验式的标称适用上限）
%     allowColdBrineExtrapolation 默认 false；仅在复现论文图3.3的
%       Tsurf=-35 degC 场景时设为 true，允许F1冷段多项式外推至-35 degC
%     interfacePermittivityMode 默认 "midpoint"（历史口径）；"boundary"
%       在空气/雪、雪/冰、冰/水物理界面用界面两侧温盐求介电常数。
%     radiativeTransferMode 默认 "burke"；"upperReflection" 传播上层
%       下行源项，用于检查层数收敛性，不作为论文复现的默认口径。
%
%   [TB, DIAGNOSTICS] 还返回层温、盐度和卤水体积分数的摘要。实现没有把
%   论文验证所得约 1.6--1.8 K 偏差硬编码为大气或经验校正项。
%
%   说明：公开论文没有给出原作者程序及完整离散代码。本函数严格使用论文
%   给出的廓线和介电参数，并用 Burke (1979) 的一阶非相干递推补全层间传输；
%   因而是可审计复现实现，而不是原作者代码的逐行复刻。
%
%   文献与定位：
%   [DS17] Zhou et al., Data Synergy..., Sect. 2.1, Eq. (1), Table 1.
%   [ZR17] Zhou et al., Improving L-band radiation model..., Sect. 3,
%          Eqs. (1)--(3), Fig. 1; Sect. 4.1.
%   [TH]   周璐博士论文, Sect. 3.2.2, Eqs. (3-1)--(3-13), Figs. 3.2--3.3.
%   [M13]  Maaß et al., Snow thickness retrieval..., Sect. 2.1,
%          Eqs. (1)--(13).
%   [K10]  Kaleschke et al. (2010), Sects. 2, 2.1, 2.2.1,
%          Eqs. (7)--(8), Tables 1--2:
%          https://doi.org/10.5194/tc-4-583-2010
%   [Z18]  Zhou et al. (2018), Sect. 2.3 and Appendix A1--A3:
%          https://doi.org/10.5194/tc-12-993-2018
%   [B79]  Burke et al. (1979), stratified incoherent emission model:
%          https://doi.org/10.1029/JC084iC01p00287
%   [CW74] Cox and Weeks (1974), ice salinity versus thickness:
%          https://doi.org/10.3189/S0022143000023418
%   [CW83] Cox and Weeks (1983), Eqs. (4)--(7), Tables I--II:
%          https://doi.org/10.1017/S0022143000008364
%   [S59]  Schwarzacher (1959), observed MYI salinity profile:
%          https://doi.org/10.1029/JZ064i012p02357
%   [T84]  Tiuri et al. (1984), snow complex permittivity:
%          https://doi.org/10.1109/JOE.1984.1145645
%   [V78]  Vant et al. (1978), sea-ice complex permittivity:
%          https://doi.org/10.1063/1.325018
%   [KS77] Klein and Swift (1977), seawater Debye model:
%          https://doi.org/10.1109/TAP.1977.1141539
%   [U64]  Untersteiner (1964), sea-ice heat conduction:
%          https://doi.org/10.1029/JZ069i022p04755
%   [YR96] Yu and Rothrock (1996), adopted thermal conductivities:
%          https://doi.org/10.1029/96JC02242

    required = {'surfaceTempC', 'waterTempC', 'waterSalinity'};
    if ~isstruct(model) || ~all(isfield(model, required))
        error('tb_forward:Parameters', ...
            'model 需包含 surfaceTempC、waterTempC、waterSalinity。');
    end
    if ~isequal(size(hi), size(hs)) || ...
            ~isnumeric(hi) || ~isnumeric(hs) || ...
            any(~isfinite(hi(:)) | hi(:) < 0) || ...
            any(~isfinite(hs(:)) | hs(:) < 0) || ...
            ~isscalar(isMYI) || ~ismember(double(isMYI), [0 1])
        error('tb_forward:Inputs', ...
            'hi、hs 需为同形状非负有限数值；isMYI 需为 0 或 1。');
    end

    model = apply_defaults(model);
    validate_model(model);

    tb = nan(size(hi));
    if nargout > 1
        diagnostics = repmat(struct( ...
            'modelName', "Zhou2017_multilayer_Burke", ...
            'nIceLayers', model.nIceLayers, ...
            'nSnowLayers', model.nSnowLayers, ...
            'bulkIceSalinityPpt', NaN, ...
            'minIceTemperatureC', NaN, ...
            'maxIceTemperatureC', NaN, ...
            'maxBrineVolumePermille', NaN, ...
            'brineValidityExceeded', false, ...
            'emissivityMean', NaN, ...
            'profile', struct(), ...
            'boundary', struct(), ...
            'interfacePermittivityMode', "midpoint", ...
            'radiativeTransferMode', "burke"), size(hi));
    else
        diagnostics = struct([]);
    end

    for index = 1:numel(hi)
        [tb(index), detail] = forward_scalar( ...
            double(hi(index)), double(hs(index)), logical(isMYI), model);
        if nargout > 1
            diagnostics(index) = detail;
        end
    end
end

function model = apply_defaults(model)
    defaults = struct( ...
        'nIceLayers', 20, ...
        'nSnowLayers', 1, ...
        'frequencyHz', 1.4e9, ...
        'incidenceAnglesDeg', 0:1:40, ...
        'includeCosmic', true, ...
        'cosmicTBK', 2.7, ...
        'atmosphereTBK', 0, ...
        'tbOffsetK', 0, ...
        'snowDensityKgM3', 320, ...
        'snowConductivity', 0.31, ...
        'fyiMinSalinityPpt', 1.5, ...
        'maxBrinePermille', 70, ...
        'allowColdBrineExtrapolation', false, ...
        'interfacePermittivityMode', "midpoint", ...
        'radiativeTransferMode', "burke");
    % 参数来源：N_ice=20、N_snow=1 见 [ZR17] Sect. 3 与 [TH] Sect. 3.2.2；
    % f=1.4 GHz 见 [DS17] Sect. 2.1；0--40 deg 双极化平均见 [ZR17] Sect. 2
    % 和 [Z18] Sect. 2.3。1 deg 步长是本代码的数值离散，并非论文指定。
    % TB_cosm=2.7 K 见 [M13] Sect. 2.1 与 [Z18] Appendix A3 Eq. (A1)。
    % rho_snow=320 kg/m^3 见 [DS17] Sect. 2.1、[ZR17] Sect. 3；
    % k_snow=0.31 W/(m K) 见 [M13] Eq. (9) 前文、[ZR17] Eq. (2)。
    % FYI 盐度下限 1.5 ppt 见 [ZR17] Sect. 3；Vb<70 permille 的适用范围
    % 见 [M13] Sect. 2.1 与 [K10] Sect. 2.2.1 Eq. (7)。
    % atmosphereTBK=0 与 tbOffsetK=0 是透明的实现默认值；[Z18] Appendix A3
    % 给出 DeltaTB_atm 的范围但未给每个场景的确定值，故不擅自校正。
    names = fieldnames(defaults);
    for index = 1:numel(names)
        if ~isfield(model, names{index}) || isempty(model.(names{index}))
            model.(names{index}) = defaults.(names{index});
        end
    end
end

function validate_model(model)
    scalarFinite = {'surfaceTempC','waterTempC','waterSalinity', ...
        'frequencyHz','cosmicTBK','atmosphereTBK','tbOffsetK', ...
        'snowDensityKgM3','snowConductivity','fyiMinSalinityPpt'};
    scalarFinite{end+1} = 'maxBrinePermille';
    for index = 1:numel(scalarFinite)
        value = model.(scalarFinite{index});
        if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value)
            error('tb_forward:Parameters', ...
                'model.%s 必须是有限标量。', scalarFinite{index});
        end
    end
    if model.surfaceTempC >= 0 || model.waterTempC >= 0 || ...
            model.waterTempC <= model.surfaceTempC || ...
            model.waterSalinity <= 0 || model.frequencyHz <= 0 || ...
            model.snowDensityKgM3 <= 0 || model.snowConductivity <= 0 || ...
            model.fyiMinSalinityPpt < 0 || model.maxBrinePermille <= 0
        error('tb_forward:Parameters', ...
            '温度、盐度、频率、雪密度或热导率超出物理范围。');
    end
    if ~isscalar(model.nIceLayers) || model.nIceLayers < 1 || ...
            model.nIceLayers ~= round(model.nIceLayers) || ...
            ~isscalar(model.nSnowLayers) || model.nSnowLayers < 1 || ...
            model.nSnowLayers ~= round(model.nSnowLayers)
        error('tb_forward:Parameters', '冰雪层数必须为正整数。');
    end
    if ~isnumeric(model.incidenceAnglesDeg) || ...
            isempty(model.incidenceAnglesDeg) || ...
            any(~isfinite(model.incidenceAnglesDeg(:))) || ...
            any(model.incidenceAnglesDeg(:) < 0 | ...
                model.incidenceAnglesDeg(:) > 40)
        error('tb_forward:Parameters', ...
            'incidenceAnglesDeg 必须位于 0--40 度。');
    end
    if ~(islogical(model.includeCosmic) || ...
            (isnumeric(model.includeCosmic) && isscalar(model.includeCosmic) && ...
            ismember(double(model.includeCosmic), [0 1])))
        error('tb_forward:Parameters', 'includeCosmic 必须为逻辑标量。');
    end
    if ~(islogical(model.allowColdBrineExtrapolation) || ...
            (isnumeric(model.allowColdBrineExtrapolation) && ...
            isscalar(model.allowColdBrineExtrapolation) && ...
            ismember(double(model.allowColdBrineExtrapolation), [0 1])))
        error('tb_forward:Parameters', ...
            'allowColdBrineExtrapolation 必须为逻辑标量。');
    end
    if ~isscalar(string(model.interfacePermittivityMode)) || ...
            ~ismember(string(model.interfacePermittivityMode), ...
            ["midpoint", "boundary"])
        error('tb_forward:Parameters', ...
            'interfacePermittivityMode 只能为 midpoint 或 boundary。');
    end
    if ~isscalar(string(model.radiativeTransferMode)) || ...
            ~ismember(string(model.radiativeTransferMode), ...
            ["burke", "upperReflection"])
        error('tb_forward:Parameters', ...
            'radiativeTransferMode 只能为 burke 或 upperReflection。');
    end
end

function [tbMean, detail] = forward_scalar(hi, hs, isMYI, model)
    if hi == 0 && hs > 0
        error('tb_forward:SnowWithoutIce', ...
            'hi=0 时不能存在非零雪厚。');
    end

    epsWater = seawater_permittivity( ...
        model.frequencyHz, model.waterTempC, model.waterSalinity);
    % 海水复介电常数：Klein-Swift 经验 Debye 模型，[KS77]；采用关系的
    % 文献链见 [M13] Sect. 2.1、[TH] Sect. 3.2.2、[Z18] Appendix A3。
    waterTempK = model.waterTempC + 273.15;

    layerThickness = [];
    layerTemperatureK = [];
    layerPermittivity = [];
    layerKind = [];
    iceTemperatureC = [];
    iceSalinityPpt = [];
    brinePermille = [];
    epsSnow = [];
    profile = struct();
    boundary = struct();

    if hi > 0
        [iceTemperatureC, snowTemperatureC, iceSalinityPpt, profile] = ...
            vertical_profiles(hi, hs, isMYI, model);
        brinePermille = brine_volume_permille( ...
            iceTemperatureC, iceSalinityPpt, ...
            model.allowColdBrineExtrapolation);
        epsIce = ice_permittivity( ...
            min(brinePermille, model.maxBrinePermille), isMYI);
        topBrine = brine_volume_permille( ...
            profile.interfaceTemperatureC, profile.iceTopSalinityPpt, ...
            model.allowColdBrineExtrapolation);
        bottomBrine = brine_volume_permille( ...
            model.waterTempC, profile.iceBottomSalinityPpt, ...
            model.allowColdBrineExtrapolation);
        boundary.iceTop = ice_permittivity( ...
            min(topBrine, model.maxBrinePermille), isMYI);
        boundary.iceBottom = ice_permittivity( ...
            min(bottomBrine, model.maxBrinePermille), isMYI);
        boundary.iceWaterLower = epsWater;
        boundary.iceTopBrineRawPermille = topBrine;
        boundary.iceBottomBrineRawPermille = bottomBrine;
        % Vb 及 eps_ice 的公式分别见 [CW83] Eqs. (5)--(7), Tables I--II
        % 和 [V78]；复现链见 [TH] Eqs. (3-9)--(3-13)。Vb 按论文直接使用
        % rho_i*S_i/F1(T)，不再用 F2 引入 Va/V=0 的额外密度闭合。
        % min(...,70) 是对 [V78]/[K10] 介电经验式有效域的保护。

        if hs > 0
            epsSnow = snow_permittivity( ...
                model.frequencyHz, snowTemperatureC, model.snowDensityKgM3);
            layerThickness = [repmat(hs / model.nSnowLayers, 1, model.nSnowLayers), ...
                repmat(hi / model.nIceLayers, 1, model.nIceLayers)];
            layerTemperatureK = [snowTemperatureC + 273.15, ...
                iceTemperatureC + 273.15];
            layerPermittivity = [epsSnow, epsIce];
            layerKind = [ones(1, model.nSnowLayers), ...
                2 * ones(1, model.nIceLayers)];
            boundary.snowTop = snow_permittivity( ...
                model.frequencyHz, model.surfaceTempC, ...
                model.snowDensityKgM3);
            boundary.snowBottom = snow_permittivity( ...
                model.frequencyHz, profile.interfaceTemperatureC, ...
                model.snowDensityKgM3);
            boundary.top = boundary.snowTop;
        else
            layerThickness = repmat(hi / model.nIceLayers, 1, model.nIceLayers);
            layerTemperatureK = iceTemperatureC + 273.15;
            layerPermittivity = epsIce;
            layerKind = 2 * ones(1, model.nIceLayers);
            boundary.top = boundary.iceTop;
        end
        profile.iceLayerBrineRawPermille = brinePermille;
        profile.iceLayerBrineUsedPermille = ...
            min(brinePermille, model.maxBrinePermille);
        profile.iceLayerPermittivity = epsIce;
        profile.snowLayerPermittivity = epsSnow;
    end

    angles = model.incidenceAnglesDeg(:).';
    samples = nan(2, numel(angles));
    emissivity = nan(2, numel(angles));
    for polarization = 1:2
        for angleIndex = 1:numel(angles)
            if string(model.radiativeTransferMode) == "upperReflection"
                [surfaceTB, effectiveEmissivity] = ...
                    radiative_stack(angles(angleIndex), polarization, ...
                    model.frequencyHz, layerThickness, layerTemperatureK, ...
                    layerPermittivity, layerKind, waterTempK, epsWater, ...
                    boundary, model.interfacePermittivityMode);
            else
                [surfaceTB, effectiveEmissivity] = ...
                    burke_stack(angles(angleIndex), polarization, ...
                    model.frequencyHz, layerThickness, layerTemperatureK, ...
                    layerPermittivity, layerKind, waterTempK, epsWater, ...
                    boundary, model.interfacePermittivityMode);
            end
            if model.includeCosmic
                surfaceTB = surfaceTB + ...
                    (1 - effectiveEmissivity) * model.cosmicTBK;
            end
            % 纯冰覆盖(c=1)时，以上为 [Z18] Appendix A3 Eq. (A1) 中的
            % TB_ice+(1-e_ice)TB_cosm；开阔水/lead 已按 [DS17] Sect. 2.1
            % 的理论场景限制排除，故此处不计算面积混合项。
            samples(polarization, angleIndex) = surfaceTB;
            emissivity(polarization, angleIndex) = effectiveEmissivity;
        end
    end
    tbMean = mean(samples, 'all') + model.atmosphereTBK + model.tbOffsetK;
    % 对 0--40 deg 及 H/V 作算术平均：见 [ZR17] Sect. 2、[Z18] Sect. 2.3。
    % atmosphereTBK 对应 [Z18] Eq. (A1) 的 DeltaTB_atm；tbOffsetK 仅预留
    % 显式标定，默认均为 0，未把论文验证偏差当作物理修正。

    if isempty(iceSalinityPpt)
        bulkSalinity = NaN;
        minIceTemperature = NaN;
        maxIceTemperature = NaN;
        maxBrine = NaN;
        brineValidityExceeded = false;
    else
        bulkSalinity = mean(iceSalinityPpt);
        minIceTemperature = min(iceTemperatureC);
        maxIceTemperature = max(iceTemperatureC);
        maxBrine = max(brinePermille);
        brineValidityExceeded = any(brinePermille >= model.maxBrinePermille);
    end
    detail = struct( ...
        'modelName', "Zhou2017_multilayer_Burke", ...
        'nIceLayers', model.nIceLayers, ...
        'nSnowLayers', model.nSnowLayers, ...
        'bulkIceSalinityPpt', bulkSalinity, ...
        'minIceTemperatureC', minIceTemperature, ...
        'maxIceTemperatureC', maxIceTemperature, ...
        'maxBrineVolumePermille', maxBrine, ...
        'brineValidityExceeded', brineValidityExceeded, ...
        'emissivityMean', mean(emissivity, 'all'), ...
        'profile', profile, ...
        'boundary', boundary, ...
        'interfacePermittivityMode', ...
        string(model.interfacePermittivityMode), ...
        'radiativeTransferMode', ...
        string(model.radiativeTransferMode));
end

%盐度廓线和温度廓线的计算
function [iceTemperatureC, snowTemperatureC, salinityPpt, profile] = ...
        vertical_profiles(hi, hs, isMYI, model)
    nIce = model.nIceLayers;
    nSnow = model.nSnowLayers;

    if isMYI
        z = ((1:nIce) - 0.5) / nIce;
        salinityPpt = 0.5 * 3.2 .* ...
            (1 - cos(pi .* z .^ (0.407 ./ (z + 0.573))));
        topSalinityPpt = 0;
        bottomSalinityPpt = 3.2;
        % MYI: z_k=(k-1/2)/N, Smax=3.2 ppt, a=0.407, b=0.573；
        % 见 [ZR17] Sect. 3 Eq. (3)、[TH] Sect. 3.2.2 Eq. (3-6)，
        % 原始观测来源为 [S59]。
    else
        bulkSalinity = 6.08 * exp(-5.81 * hi) + ...
            7.409 * exp(-0.5228 * hi);
        bulkSalinity = max(model.fyiMinSalinityPpt, bulkSalinity);
        salinityPpt = repmat(bulkSalinity, 1, nIce);
        topSalinityPpt = bulkSalinity;
        bottomSalinityPpt = bulkSalinity;
        % FYI: S=6.08 exp(-5.81 h_i)+7.409 exp(-0.5228 h_i)，各层同盐度；
        % 见 [TH] Sect. 3.2.2 Eq. (3-5)，文献源为 [CW74]；1.5 ppt 下限
        % 见 [ZR17] Sect. 3 与 [TH] Eq. (3-5) 后说明。
    end

    meanColumnTempC = 0.5 * (model.surfaceTempC + model.waterTempC);
    kIce = 2.034 + 0.13 * mean(salinityPpt) / meanColumnTempC;
    % k_ice=2.034+0.13*S_ice/(T_ice-273.15)，此处分母直接用摄氏温度；
    % 见 [M13] Sect. 2.1 Eq. (9)、[ZR17] Sect. 3 Eq. (1)、[TH] Eq. (3-1)。
    % 采用列平均温度替代未知 bulk T_ice 的处理见 [M13] Eq. (9) 后说明；
    % 物理依据追溯至 [U64]，参数采用见 [YR96]。
    if kIce <= 0
        error('tb_forward:ThermalConductivity', ...
            '计算所得海冰热导率非正。');
    end

    if hs > 0
        totalResistance = hs / model.snowConductivity + hi / kIce;
        heatFlux = (model.waterTempC - model.surfaceTempC) / totalResistance;
        interfaceTempC = model.surfaceTempC + ...
            heatFlux * hs / model.snowConductivity;
        snowDepth = (((1:nSnow) - 0.5) / nSnow) * hs;
        snowTemperatureC = model.surfaceTempC + ...
            heatFlux .* snowDepth / model.snowConductivity;
    else
        interfaceTempC = model.surfaceTempC;
        snowTemperatureC = [];
        snowDepth = [];
        heatFlux = (model.waterTempC - model.surfaceTempC) / (hi / kIce);
    end
    iceDepth = (((1:nIce) - 0.5) / nIce) * hi;
    iceTemperatureC = interfaceTempC + heatFlux .* iceDepth / kIce;
    profile = struct( ...
        'surfaceTemperatureC', model.surfaceTempC, ...
        'interfaceTemperatureC', interfaceTempC, ...
        'snowBottomTemperatureC', ...
        model.surfaceTempC + heatFlux * hs / model.snowConductivity, ...
        'iceBottomTemperatureC', ...
        interfaceTempC + heatFlux * hi / kIce, ...
        'waterTemperatureC', model.waterTempC, ...
        'snowLayerMidDepthM', snowDepth, ...
        'iceLayerMidDepthM', iceDepth, ...
        'snowLayerTemperatureC', snowTemperatureC, ...
        'iceLayerTemperatureC', iceTemperatureC, ...
        'iceLayerSalinityPpt', salinityPpt, ...
        'iceTopSalinityPpt', topSalinityPpt, ...
        'iceBottomSalinityPpt', bottomSalinityPpt, ...
        'iceConductivityWmK', kIce, ...
        'snowConductivityWmK', model.snowConductivity, ...
        'heatFluxWm2', heatFlux, ...
        'bulkIceTemperatureApproxC', meanColumnTempC);
    % 稳态一维导热：R=h_s/k_s+h_i/k_i, q=(T_w-T_surf)/R，冰雪界面
    % 温度连续，各层取中点温度。见 [M13] Sect. 2.1 Eqs. (4)--(13)、
    % [ZR17] Sect. 3 Eqs. (1)--(2)、[TH] Sect. 3.2.2 Eqs. (3-1)--(3-4)。
end

function epsSnow = snow_permittivity(frequencyHz, temperatureC, densityKgM3)
    relativeDensity = densityKgM3 / 1000;
    realPart = 1 + 0.7 * relativeDensity + 0.7 * relativeDensity^2;
    lossPart = 1.59e6 * ...
        (0.52 * relativeDensity + 0.62 * relativeDensity^2) .* ...
        (frequencyHz^-1 + 1.23e-14 * sqrt(frequencyHz)) .* ...
        exp(0.036 .* temperatureC);
    epsSnow = realPart + 1i .* lossPart;
    % 干雪复介电常数拟合：见 [TH] Sect. 3.2.2 Eq. (3-8) 与 [Z18]
    % Appendix A3；原始实验文献为 [T84]。rho 为相对水密度，f 为 Hz，T 为 degC。
end

% 计算海冰介电常数
function epsIce = ice_permittivity(brinePermille, isMYI)
    if isMYI
        coefficients = [3.10, 0.0084, 0.003, 0.00435];
    else
        coefficients = [3.10, 0.0084, 0.037, 0.00445];
    end
    epsIce = coefficients(1) + coefficients(2) .* brinePermille + ...
        1i .* (coefficients(3) + coefficients(4) .* brinePermille);
    % eps_ice=a1+a2*Vb+i(a3+a4*Vb)：[V78]；1.4 GHz FYI/MYI 系数取自
    % [K10] Sect. 2.2.1 Eq. (7), Table 1；[TH] Eq. (3-9)。Vb 单位为 permille。
end

% 卤水体积分数计算，用于计算介电常数
function vbPermille = brine_volume_permille( ...
        temperatureC, salinityPpt, allowColdExtrapolation)
    if nargin < 3
        allowColdExtrapolation = false;
    end
    minimumTemperatureC = -30.000001;
    if allowColdExtrapolation
        minimumTemperatureC = -35.000001;
    end
    if any(temperatureC >= 0) || any(temperatureC < minimumTemperatureC)
        error('tb_forward:BrineTemperature', ...
            ['卤水体积分数公式要求冰层温度位于 -30--0 摄氏度；', ...
            '论文图3.3的显式外推选项最多允许至 -35 摄氏度。']);
    end
    f1 = nan(size(temperatureC));
    warm = temperatureC > -2;
    middle = temperatureC <= -2 & temperatureC > -22.9;
    cold = temperatureC <= -22.9;
    f1(warm) = polyval([0.21454, 0.58402, -18.407, -0.041221], ...
        temperatureC(warm));
    % -2<T<0 degC 系数：Lepparanta-Manninen 参数化，见 [K10]
    % Sect. 2.2.1 Eq. (8), Table 2。
    f1(middle) = polyval([-0.01074, -0.6397, -22.45, -4.732], ...
        temperatureC(middle));
    f1(cold) = polyval([0.7160, 55.27, 1309, 9899], ...
        temperatureC(cold));
    % -30<=T<=-2 degC 的分段 F1 系数：见 [CW83] Eqs. (5)--(7),
    % Table II；F1=rho_b*S_b*(1+k)，-22.9 degC 是固体盐相变分段点。
    % Cox--Weeks表值到-30 degC。周璐论文图3.3明确使用Tsurf=-35 degC，
    % 因而仅在allowColdBrineExtrapolation=true时将该冷段多项式外推至
    % -35 degC；默认计算仍严格限制在原表范围内。

    pureIceDensity = 0.917 - 1.403e-4 .* temperatureC;
    vbFraction = pureIceDensity .* salinityPpt ./ f1;
    % rho_i=0.917-1.403e-4*T，Vb=rho_i*S_i/F1(T)，随后乘 1000 转为
    % permille；见 [CW83] Eqs. (5)--(7), Tables I--II 与 [TH]
    % Eqs. (3-10)--(3-13)。
    %
    % 旧实现还计算 F2，并等价地令气体相对体积 Va/V=0 后，用 [CW83]
    % Eq. (15) 消去整体海冰密度 rho。该额外闭合不属于论文直接步骤，
    % 故注释保留而不再参与默认计算：
    % f2 = nan(size(temperatureC));
    % f2(warm) = polyval([1.3603e-4, 1.2291e-4, -0.016111, 0.090312], ...
    %     temperatureC(warm));
    % f2(middle) = polyval([-8.801e-6, -5.330e-4, -0.01763, 0.08903], ...
    %     temperatureC(middle));
    % f2(cold) = polyval([5.819e-4, 0.04518, 1.089, 8.547], ...
    %     temperatureC(cold));
    % denominator = f1 - pureIceDensity .* salinityPpt .* f2;
    % vbFraction = pureIceDensity .* salinityPpt ./ denominator;
    if any(~isfinite(vbFraction) | vbFraction < 0)
        error('tb_forward:BrineVolume', ...
            '卤水体积分数计算失败；请检查温度和盐度。');
    end
    vbPermille = 1000 .* vbFraction;
end

function epsWater = seawater_permittivity(frequencyHz, temperatureC, salinityPpt)
    salinity2 = salinityPpt^2;
    salinity3 = salinityPpt^3;
    temperature2 = temperatureC^2;
    temperature3 = temperatureC^3;
    staticPermittivity = 87.134 - 0.1949 * temperatureC - ...
        0.01276 * temperature2 + 2.491e-4 * temperature3;
    staticScale = 1 + 1.613e-5 * salinityPpt * temperatureC - ...
        3.656e-3 * salinityPpt + 3.21e-5 * salinity2 - ...
        4.232e-7 * salinity3;
    staticPermittivity = staticPermittivity * staticScale;

    relaxationTime = 1.768e-11 - 6.086e-13 * temperatureC + ...
        1.104e-14 * temperature2 - 8.111e-17 * temperature3;
    relaxationScale = 1 + 2.282e-5 * salinityPpt * temperatureC - ...
        7.638e-4 * salinityPpt - 7.760e-6 * salinity2 + ...
        1.105e-8 * salinity3;
    relaxationTime = relaxationTime * relaxationScale;

    conductivity = salinityPpt * (0.182521 - ...
        1.46192e-3 * salinityPpt + 2.09324e-5 * salinity2 - ...
        1.28205e-7 * salinity3);
    delta = 25 - temperatureC;
    beta = 2.033e-2 + 1.266e-4 * delta + ...
        2.464e-6 * delta^2 - salinityPpt * ...
        (1.849e-5 - 2.551e-7 * delta + 2.551e-8 * delta^2);
    conductivity = conductivity * exp(-beta * delta);

    omega = 2 * pi * frequencyHz;
    epsilon0 = 8.854e-12;
    raw = 4.9 + (staticPermittivity - 4.9) / ...
        (1 + 1i * omega * relaxationTime) - ...
        1i * conductivity / (omega * epsilon0);
    epsWater = real(raw) + 1i * (-imag(raw));
    % 单 Debye 海水模型：eps=eps_inf+(eps_s-eps_inf)/(1+i*omega*tau)
    % -i*sigma/(omega*eps0)，eps_inf=4.9；eps_s、tau、sigma 的温盐
    % 多项式均来自 [KS77]。采用该关系见 [M13] Sect. 2.1、[Z18] Appendix A3。
end
