function result = audit_layer_contributions_scenario_viii(outputDirectory)
%AUDIT_LAYER_CONTRIBUTIONS_SCENARIO_VIII Decompose scenario VIII brightness.
%   Uses the current model convention: physical interfaces retain one
%   Fresnel reflection, while numerical interfaces inside continuous ice
%   have zero reflection. Five-layer blocks are exact affine compositions,
%   not averaged replacement layers.

    if nargin < 1 || strlength(string(outputDirectory)) == 0
        outputDirectory = fullfile(fileparts(mfilename('fullpath')), ...
            'layer_contribution_audit');
    end
    if ~isfolder(outputDirectory), mkdir(outputDirectory); end

    matlabRoot = fileparts(fileparts(mfilename('fullpath')));
    addpath(matlabRoot);
    data = readtable(fullfile(matlabRoot, 'processed', ...
        'tb_forward_table1_validation.csv'));
    row = data(string(data.scenario) == "VIII", :);
    assert(height(row) == 1, 'Expected one scenario VIII row.');

    model = struct('surfaceTempC', -30, 'waterTempC', -1.8, ...
        'waterSalinity', 33);
    [tbDirect, detail] = tb_forward(row.hi, row.hs, row.isMYI, model);

    epsSnow = detail.profile.snowLayerPermittivity(:);
    epsIce = detail.profile.iceLayerPermittivity(:);
    tempSnowK = detail.profile.snowLayerTemperatureC(:) + 273.15;
    tempIceK = detail.profile.iceLayerTemperatureC(:) + 273.15;
    nSnow = numel(epsSnow);
    nIce = numel(epsIce);
    assert(nSnow == 1 && nIce == 20, ...
        'This audit expects the default 1 snow + 20 ice layers.');

    epsLayer = [epsSnow; epsIce];
    tempLayerK = [tempSnowK; tempIceK];
    kind = [ones(nSnow,1); 2 * ones(nIce,1)];
    thicknessM = [row.hs / nSnow; repmat(row.hi / nIce, nIce, 1)];
    nLayer = numel(epsLayer);
    waterTempK = model.waterTempC + 273.15;
    epsWater = detail.boundary.iceWaterLower;
    anglesDeg = (0:40).';
    frequencyHz = 1.4e9;
    cosmicTBK = 2.7;

    nSamples = 2 * numel(anglesDeg);
    surfaceContribution = nan(nSamples, nLayer + 1); % water, then layers
    localSource = nan(nSamples, nLayer);
    oneWayTau = nan(nSamples, nLayer);
    escapeFactor = nan(nSamples, nLayer + 1);
    cosmicContribution = nan(nSamples, 1);
    blockTB = nan(nSamples, 1);
    sample = 0;

    for polarization = 1:2
        for angleIndex = 1:numel(anglesDeg)
            sample = sample + 1;
            sinTheta0 = sind(anglesDeg(angleIndex));
            k0 = 2 * pi * frequencyHz / 299792458;
            lowerPermittivity = epsWater;
            lowerKind = 3;
            layerA = nan(nLayer,1);
            layerS = nan(nLayer,1);
            layerE = nan(nLayer,1);

            for layer = nLayer:-1:1
                epsCurrent = epsLayer(layer);
                if kind(layer) == lowerKind
                    reflectionBottom = 0;
                else
                    reflectionBottom = fresnel_power_local( ...
                        epsCurrent, lowerPermittivity, sinTheta0, polarization);
                end
                verticalIndex = sqrt(epsCurrent - sinTheta0^2);
                tau = exp(-2 * k0 * abs(imag(verticalIndex)) * ...
                    thicknessM(layer));
                sourceFactor = (1 - tau) * (1 + reflectionBottom * tau);
                layerA(layer) = (1 - reflectionBottom) * tau;
                layerS(layer) = tempLayerK(layer) * sourceFactor;
                layerE(layer) = sourceFactor;
                oneWayTau(sample, layer) = tau;
                localSource(sample, layer) = layerS(layer);
                lowerPermittivity = epsCurrent;
                lowerKind = kind(layer);
            end

            reflectionTop = fresnel_power_local( ...
                1, epsLayer(1), sinTheta0, polarization);
            topTransmission = 1 - reflectionTop;

            contributions = zeros(nLayer + 1, 1);
            contributions(1) = waterTempK;
            emissivityUnit = 1;
            for layer = nLayer:-1:1
                contributions = layerA(layer) .* contributions;
                contributions(layer + 1) = layerS(layer);
                emissivityUnit = layerE(layer) + ...
                    layerA(layer) * emissivityUnit;
            end
            contributions = topTransmission .* contributions;
            surfaceContribution(sample,:) = contributions.';
            escapeFactor(sample,1) = contributions(1) / waterTempK;
            for layer = 1:nLayer
                escapeFactor(sample,layer+1) = ...
                    contributions(layer+1) / layerS(layer);
            end
            effectiveEmissivity = topTransmission * emissivityUnit;
            cosmicContribution(sample) = ...
                (1 - effectiveEmissivity) * cosmicTBK;

            % Exact grouping of the 20 ice layers into four affine blocks.
            upward = waterTempK;
            upwardUnit = 1;
            iceStart = nSnow + 1;
            for block = 4:-1:1
                firstLayer = iceStart + (block - 1) * 5;
                lastLayer = firstLayer + 4;
                [groupS, groupA] = compose_affine( ...
                    layerS(firstLayer:lastLayer), ...
                    layerA(firstLayer:lastLayer));
                [groupE, groupAE] = compose_affine( ...
                    layerE(firstLayer:lastLayer), ...
                    layerA(firstLayer:lastLayer));
                upward = groupS + groupA * upward;
                upwardUnit = groupE + groupAE * upwardUnit;
            end
            for layer = nSnow:-1:1
                upward = layerS(layer) + layerA(layer) * upward;
                upwardUnit = layerE(layer) + layerA(layer) * upwardUnit;
            end
            blockTB(sample) = topTransmission * upward + ...
                (1 - topTransmission * upwardUnit) * cosmicTBK;
        end
    end

    meanContribution = mean(surfaceContribution, 1).';
    meanCosmic = mean(cosmicContribution);
    tbFromComponents = sum(meanContribution) + meanCosmic;
    tbFromExactBlocks = mean(blockTB);
    assert(abs(tbFromComponents - tbDirect) < 1e-9, ...
        'Layer components do not close to tb_forward.');
    assert(abs(tbFromExactBlocks - tbDirect) < 1e-9, ...
        'Exact five-layer blocks do not close to layer-by-layer recurrence.');

    depthTopM = [0; cumsum(thicknessM(1:end-1))];
    depthBottomM = cumsum(thicknessM);
    depthCenterM = 0.5 * (depthTopM + depthBottomM);
    material = [repmat("snow",nSnow,1); repmat("ice",nIce,1)];
    layerInMaterial = [(1:nSnow).'; (1:nIce).'];
    brinePermille = [nan(nSnow,1); ...
        detail.profile.iceLayerBrineUsedPermille(:)];

    layerTable = table(material, layerInMaterial, depthTopM, depthBottomM, ...
        depthCenterM, thicknessM, tempLayerK, brinePermille, ...
        real(epsLayer), imag(epsLayer), mean(oneWayTau,1).', ...
        mean(localSource,1).', mean(escapeFactor(:,2:end),1).', ...
        meanContribution(2:end), ...
        100 .* meanContribution(2:end) ./ tbDirect, ...
        'VariableNames', {'material','layerIndex','depthTopM','depthBottomM', ...
        'depthCenterM','thicknessM','temperatureK','brineUsedPermille', ...
        'epsilonReal','epsilonImag','meanOneWayPowerTransmission', ...
        'meanLocalSourceK','meanEscapeFactorToAir', ...
        'surfaceContributionK','fractionOfTotalTBPercent'});

    extraTable = table(["water";"cosmic"], ...
        [meanContribution(1);meanCosmic], ...
        100 .* [meanContribution(1);meanCosmic] ./ tbDirect, ...
        'VariableNames', {'component','surfaceContributionK', ...
        'fractionOfTotalTBPercent'});

    blockName = ["snow";"ice_01_05";"ice_06_10"; ...
        "ice_11_15";"ice_16_20";"water";"cosmic"];
    blockDepthTopM = [0; row.hs; row.hs + 1.25; row.hs + 2.50; ...
        row.hs + 3.75; row.hs + row.hi; NaN];
    blockDepthBottomM = [row.hs; row.hs + 1.25; row.hs + 2.50; ...
        row.hs + 3.75; row.hs + row.hi; NaN; NaN];
    blockContributionK = [sum(meanContribution(2:1+nSnow)); ...
        sum(meanContribution(2+nSnow:6+nSnow)); ...
        sum(meanContribution(7+nSnow:11+nSnow)); ...
        sum(meanContribution(12+nSnow:16+nSnow)); ...
        sum(meanContribution(17+nSnow:21+nSnow)); ...
        meanContribution(1); meanCosmic];
    blockTable = table(blockName, blockDepthTopM, blockDepthBottomM, ...
        blockContributionK, 100 .* blockContributionK ./ tbDirect, ...
        'VariableNames', {'component','depthTopM','depthBottomM', ...
        'surfaceContributionK','fractionOfTotalTBPercent'});

    layerCounts = [1;4;5;10;20;100];
    layerCountTB = nan(size(layerCounts));
    for index = 1:numel(layerCounts)
        countModel = model;
        countModel.nIceLayers = layerCounts(index);
        layerCountTB(index) = tb_forward( ...
            row.hi, row.hs, row.isMYI, countModel);
    end
    layerCountTable = table(layerCounts, layerCountTB, ...
        layerCountTB - tbDirect, ...
        'VariableNames', {'nIceLayers','tbK','differenceFrom20LayerK'});

    methodTable = table( ...
        ["20 layers, one-by-one"; ...
         "4 exact affine blocks of 5 layers"; ...
         "4 representative layers (re-discretized)"], ...
        [tbDirect;tbFromExactBlocks;layerCountTB(layerCounts == 4)], ...
        [0;tbFromExactBlocks-tbDirect; ...
         layerCountTB(layerCounts == 4)-tbDirect], ...
        'VariableNames', {'method','tbK','differenceFrom20LayerK'});

    % Local block sensitivities. These are diagnostics, not calibration:
    % temperature, epsilon real/imaginary, or Vb is perturbed only in the
    % selected block while every other profile quantity is held fixed.
    blockLabels = ["ice_01_05";"ice_06_10";"ice_11_15"; ...
        "ice_16_20";"all_ice"];
    blockIndices = {1:5,6:10,11:15,16:20,1:20};
    nSensitivity = numel(blockIndices);
    deltaTempPlus1K = nan(nSensitivity,1);
    deltaImagPlus10Pct = nan(nSensitivity,1);
    deltaImagMinus10Pct = nan(nSensitivity,1);
    deltaRealPlus1Pct = nan(nSensitivity,1);
    deltaRealMinus1Pct = nan(nSensitivity,1);
    deltaVbPlus10Pct = nan(nSensitivity,1);
    deltaVbMinus10Pct = nan(nSensitivity,1);
    for sensitivityIndex = 1:nSensitivity
        iceIndices = blockIndices{sensitivityIndex};
        layerIndices = nSnow + iceIndices;

        testTemperature = tempLayerK;
        testTemperature(layerIndices) = ...
            testTemperature(layerIndices) + 1;
        deltaTempPlus1K(sensitivityIndex) = forward_from_layers( ...
            epsLayer, testTemperature, kind, thicknessM, waterTempK, ...
            epsWater, anglesDeg, frequencyHz, cosmicTBK) - tbDirect;

        testPermittivity = epsLayer;
        testPermittivity(layerIndices) = ...
            real(testPermittivity(layerIndices)) + ...
            1i * 1.10 .* imag(testPermittivity(layerIndices));
        deltaImagPlus10Pct(sensitivityIndex) = forward_from_layers( ...
            testPermittivity, tempLayerK, kind, thicknessM, waterTempK, ...
            epsWater, anglesDeg, frequencyHz, cosmicTBK) - tbDirect;
        testPermittivity = epsLayer;
        testPermittivity(layerIndices) = ...
            real(testPermittivity(layerIndices)) + ...
            1i * 0.90 .* imag(testPermittivity(layerIndices));
        deltaImagMinus10Pct(sensitivityIndex) = forward_from_layers( ...
            testPermittivity, tempLayerK, kind, thicknessM, waterTempK, ...
            epsWater, anglesDeg, frequencyHz, cosmicTBK) - tbDirect;

        testPermittivity = epsLayer;
        testPermittivity(layerIndices) = ...
            1.01 .* real(testPermittivity(layerIndices)) + ...
            1i .* imag(testPermittivity(layerIndices));
        deltaRealPlus1Pct(sensitivityIndex) = forward_from_layers( ...
            testPermittivity, tempLayerK, kind, thicknessM, waterTempK, ...
            epsWater, anglesDeg, frequencyHz, cosmicTBK) - tbDirect;
        testPermittivity = epsLayer;
        testPermittivity(layerIndices) = ...
            0.99 .* real(testPermittivity(layerIndices)) + ...
            1i .* imag(testPermittivity(layerIndices));
        deltaRealMinus1Pct(sensitivityIndex) = forward_from_layers( ...
            testPermittivity, tempLayerK, kind, thicknessM, waterTempK, ...
            epsWater, anglesDeg, frequencyHz, cosmicTBK) - tbDirect;

        for directionIndex = 1:2
            scale = [1.10,0.90];
            vbTest = detail.profile.iceLayerBrineUsedPermille(:);
            vbTest(iceIndices) = scale(directionIndex) .* vbTest(iceIndices);
            epsIceTest = 3.10 + 0.0084 .* vbTest + ...
                1i .* (0.003 + 0.00435 .* vbTest);
            testPermittivity = [epsSnow;epsIceTest];
            deltaValue = forward_from_layers( ...
                testPermittivity, tempLayerK, kind, thicknessM, waterTempK, ...
                epsWater, anglesDeg, frequencyHz, cosmicTBK) - tbDirect;
            if directionIndex == 1
                deltaVbPlus10Pct(sensitivityIndex) = deltaValue;
            else
                deltaVbMinus10Pct(sensitivityIndex) = deltaValue;
            end
        end
    end
    paperGapK = row.tbPaper - tbDirect;
    temperatureIncreaseToCloseGapK = paperGapK ./ deltaTempPlus1K;
    sensitivityTable = table(blockLabels,deltaTempPlus1K, ...
        temperatureIncreaseToCloseGapK,deltaImagPlus10Pct, ...
        deltaImagMinus10Pct,deltaRealPlus1Pct,deltaRealMinus1Pct, ...
        deltaVbPlus10Pct,deltaVbMinus10Pct);

    writetable(layerTable, fullfile(outputDirectory, ...
        'scenario_viii_layer_contributions.csv'));
    writetable(extraTable, fullfile(outputDirectory, ...
        'scenario_viii_external_contributions.csv'));
    writetable(blockTable, fullfile(outputDirectory, ...
        'scenario_viii_five_layer_blocks.csv'));
    writetable(layerCountTable, fullfile(outputDirectory, ...
        'scenario_viii_layer_count_comparison.csv'));
    writetable(methodTable, fullfile(outputDirectory, ...
        'scenario_viii_grouping_identity.csv'));
    writetable(sensitivityTable, fullfile(outputDirectory, ...
        'scenario_viii_block_sensitivity.csv'));

    result = struct('tbDirectK',tbDirect, ...
        'tbFromComponentsK',tbFromComponents, ...
        'tbFromExactBlocksK',tbFromExactBlocks, ...
        'paperTBK',row.tbPaper, 'layerTable',layerTable, ...
        'blockTable',blockTable, 'methodTable',methodTable, ...
        'sensitivityTable',sensitivityTable);
end

function [sourceCombined, transferCombined] = compose_affine(source, transfer)
    sourceCombined = 0;
    transferCombined = 1;
    for index = numel(source):-1:1
        sourceCombined = source(index) + transfer(index) * sourceCombined;
        transferCombined = transfer(index) * transferCombined;
    end
end

function reflectivity = fresnel_power_local( ...
        epsUpper, epsLower, sinTheta0, polarization)
    qUpper = sqrt(epsUpper - sinTheta0^2);
    qLower = sqrt(epsLower - sinTheta0^2);
    if polarization == 1
        coefficient = (qUpper - qLower) / (qUpper + qLower);
    else
        coefficient = (epsLower * qUpper - epsUpper * qLower) / ...
            (epsLower * qUpper + epsUpper * qLower);
    end
    reflectivity = min(max(abs(coefficient)^2, 0), 1);
end

function tbMean = forward_from_layers(epsLayer, tempLayerK, kind, ...
        thicknessM, waterTempK, epsWater, anglesDeg, frequencyHz, cosmicTBK)
    samples = nan(2,numel(anglesDeg));
    for polarization = 1:2
        for angleIndex = 1:numel(anglesDeg)
            sinTheta0 = sind(anglesDeg(angleIndex));
            k0 = 2 * pi * frequencyHz / 299792458;
            upward = waterTempK;
            upwardUnit = 1;
            lowerPermittivity = epsWater;
            lowerKind = 3;
            for layer = numel(epsLayer):-1:1
                epsCurrent = epsLayer(layer);
                if kind(layer) == lowerKind
                    reflectionBottom = 0;
                else
                    reflectionBottom = fresnel_power_local( ...
                        epsCurrent, lowerPermittivity, ...
                        sinTheta0, polarization);
                end
                tau = exp(-2 * k0 * ...
                    abs(imag(sqrt(epsCurrent - sinTheta0^2))) * ...
                    thicknessM(layer));
                sourceFactor = (1 - tau) * ...
                    (1 + reflectionBottom * tau);
                upward = tempLayerK(layer) * sourceFactor + ...
                    (1 - reflectionBottom) * upward * tau;
                upwardUnit = sourceFactor + ...
                    (1 - reflectionBottom) * upwardUnit * tau;
                lowerPermittivity = epsCurrent;
                lowerKind = kind(layer);
            end
            reflectionTop = fresnel_power_local( ...
                1, epsLayer(1), sinTheta0, polarization);
            effectiveEmissivity = (1 - reflectionTop) * upwardUnit;
            samples(polarization,angleIndex) = ...
                (1 - reflectionTop) * upward + ...
                (1 - effectiveEmissivity) * cosmicTBK;
        end
    end
    tbMean = mean(samples,'all');
end
