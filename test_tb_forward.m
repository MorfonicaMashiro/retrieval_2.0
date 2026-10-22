function tests = test_tb_forward
%TEST_TB_FORWARD 前向亮温模型的数值与接口回归测试。
    tests = functiontests(localfunctions);
end

function testPaperTableAgreement(testCase)
    result = validate_tb_forward();
    verifyLessThan(testCase, result.scenarioRMSEK, 3.0);
    verifyLessThan(testCase, max(abs(result.scenarios.differenceK)), 4.0);
end

function testPublishedIceTypeContrast(testCase)
    model = struct('surfaceTempC', -35, 'waterTempC', -1.8, ...
        'waterSalinity', 33);
    tbFYI = tb_forward(1, 0.1, false, model);
    tbMYI = tb_forward(1, 0.1, true, model);
    % Zhou 等（2017）报告该场景 FYI 比 MYI 约高 12 K。
    verifyEqual(testCase, tbFYI - tbMYI, 12, 'AbsTol', 1.5);
end

function testThickIceLayerConvergence(testCase)
    model = struct('surfaceTempC', -30, 'waterTempC', -1.8, ...
        'waterSalinity', 33, 'nIceLayers', 20);
    tb20 = tb_forward(5, 0.4, true, model);
    model.nIceLayers = 100;
    tb100 = tb_forward(5, 0.4, true, model);
    verifyLessThan(testCase, abs(tb100 - tb20), 0.1);
end

function testVectorShapeAndFiniteOutput(testCase)
    model = struct('surfaceTempC', -30, 'waterTempC', -1.8, ...
        'waterSalinity', 33);
    hi = [0.5, 1.0; 1.5, 2.5];
    hs = [0.05, 0.03; 0.10, 0.25];
    tb = tb_forward(hi, hs, false, model);
    verifySize(testCase, tb, size(hi));
    verifyTrue(testCase, all(isfinite(tb), 'all'));
    verifyGreaterThan(testCase, min(tb, [], 'all'), 150);
    verifyLessThan(testCase, max(tb, [], 'all'), 280);
end

function testSnowWithoutIceRejected(testCase)
    model = struct('surfaceTempC', -30, 'waterTempC', -1.8, ...
        'waterSalinity', 33);
    verifyError(testCase, @() tb_forward(0, 0.1, false, model), ...
        'tb_forward:SnowWithoutIce');
end

function testColdBrineExtrapolationRequiresExplicitOptIn(testCase)
    model = struct('surfaceTempC', -35, 'waterTempC', -1.8, ...
        'waterSalinity', 33);
    verifyError(testCase, @() tb_forward(2, 0.01, false, model), ...
        'tb_forward:BrineTemperature');

    model.allowColdBrineExtrapolation = true;
    tb = tb_forward(2, 0.01, false, model);
    verifyTrue(testCase, isfinite(tb));
end

function testPhysicalBoundaryAndLayerTemperaturesAreDistinct(testCase)
    model = struct('surfaceTempC', -30, 'waterTempC', -1.8, ...
        'waterSalinity', 33);
    [tbMidpoint, detail] = tb_forward(2.5, 0.15, true, model);
    interfaceC = detail.profile.interfaceTemperatureC;
    verifyEqual(testCase, detail.profile.snowBottomTemperatureC, ...
        interfaceC, 'AbsTol', 1e-10);
    verifyGreaterThan(testCase, interfaceC, ...
        detail.profile.snowLayerTemperatureC(end));
    verifyGreaterThan(testCase, ...
        detail.profile.iceLayerTemperatureC(1), interfaceC);
    verifyNotEqual(testCase, detail.boundary.iceTop, ...
        detail.profile.iceLayerPermittivity(1));

    model.interfacePermittivityMode = "boundary";
    tbBoundary = tb_forward(2.5, 0.15, true, model);
    verifyGreaterThan(testCase, tbBoundary, tbMidpoint);
end

function testBrightnessOffsetAppliedOnce(testCase)
    model = struct('surfaceTempC', -30, 'waterTempC', -1.8, ...
        'waterSalinity', 33, 'includeCosmic', true);
    tbBase = tb_forward(1, 0.1, false, model);
    model.tbOffsetK = 1;
    verifyEqual(testCase, tb_forward(1, 0.1, false, model) - tbBase, ...
        1, 'AbsTol', 1e-10);
end

function testBrineVolumeUsesDirectF1Formula(testCase)
    model = struct('surfaceTempC', -30, 'waterTempC', -1.8, ...
        'waterSalinity', 33, 'nIceLayers', 20);
    [~, detail] = tb_forward(0.5, 0.05, false, model);
    temperatureC = detail.profile.iceLayerTemperatureC(end);
    salinityPpt = detail.profile.iceLayerSalinityPpt(end);
    f1 = polyval([-0.01074, -0.6397, -22.45, -4.732], temperatureC);
    pureIceDensity = 0.917 - 1.403e-4 * temperatureC;
    expectedPermille = 1000 * pureIceDensity * salinityPpt / f1;
    actualPermille = detail.profile.iceLayerBrineRawPermille(end);
    verifyEqual(testCase, actualPermille, expectedPermille, 'AbsTol', 1e-12);

    % 含 F2 的旧式等价于额外采用 Va/V=0 的密度闭合，不应再成为默认值。
    f2 = polyval([-8.801e-6, -5.330e-4, -0.01763, 0.08903], temperatureC);
    legacyPermille = 1000 * pureIceDensity * salinityPpt / ...
        (f1 - pureIceDensity * salinityPpt * f2);
    verifyGreaterThan(testCase, abs(actualPermille - legacyPermille), 0.1);
end

function testUniformSlabSubdivisionIsInvariant(testCase)
    frequencyHz = 1.4e9;
    thickness = 0.75;
    iceTemperatureK = 258;
    icePermittivity = 3.25 + 0.18i;
    waterTemperatureK = 271.35;
    waterPermittivity = 72 + 55i;
    boundary = struct();

    tb1 = radiative_stack(40, 1, frequencyHz, thickness, ...
        iceTemperatureK, icePermittivity, 2, ...
        waterTemperatureK, waterPermittivity, boundary, "midpoint");
    layerCount = 100;
    tb100 = radiative_stack(40, 1, frequencyHz, ...
        repmat(thickness / layerCount, 1, layerCount), ...
        repmat(iceTemperatureK, 1, layerCount), ...
        repmat(icePermittivity, 1, layerCount), ...
        2 * ones(1, layerCount), waterTemperatureK, ...
        waterPermittivity, boundary, "midpoint");
    verifyEqual(testCase, tb100, tb1, 'AbsTol', 1e-10);
end

function testUniformSlabMatchesAnalyticSolution(testCase)
    incidenceDeg = 35;
    frequencyHz = 1.4e9;
    thickness = 0.5;
    iceTemperatureK = 260;
    icePermittivity = 3.3 + 0.12i;
    waterTemperatureK = 271;
    waterPermittivity = 70 + 50i;
    boundary = struct();

    actual = radiative_stack(incidenceDeg, 1, frequencyHz, thickness, ...
        iceTemperatureK, icePermittivity, 2, ...
        waterTemperatureK, waterPermittivity, boundary, "midpoint");

    sinTheta0 = sind(incidenceDeg);
    k0 = 2 * pi * frequencyHz / 299792458;
    transmission = exp(-2 * k0 * ...
        abs(imag(sqrt(icePermittivity - sinTheta0^2))) * thickness);
    topReflection = local_fresnel_power( ...
        1, icePermittivity, sinTheta0, 1);
    bottomReflection = local_fresnel_power( ...
        icePermittivity, waterPermittivity, sinTheta0, 1);
    expected = (1 - topReflection) * ( ...
        (1 - bottomReflection) * transmission * waterTemperatureK + ...
        (1 - bottomReflection * transmission^2 - ...
        (1 - bottomReflection) * transmission) * iceTemperatureK);
    verifyEqual(testCase, actual, expected, 'AbsTol', 1e-10);
end

function testUpperReflectionModeConvergesForThinIce(testCase)
    baseModel = struct('surfaceTempC', -30, 'waterTempC', -1.8, ...
        'waterSalinity', 33, 'nIceLayers', 20);
    base20 = tb_forward(0.5, 0.05, true, baseModel);
    baseModel.nIceLayers = 100;
    base100 = tb_forward(0.5, 0.05, true, baseModel);

    model = struct('surfaceTempC', -30, 'waterTempC', -1.8, ...
        'waterSalinity', 33, 'radiativeTransferMode', "upperReflection", ...
        'nIceLayers', 20);
    tb20 = tb_forward(0.5, 0.05, true, model);
    model.nIceLayers = 100;
    [tb100, detail] = tb_forward(0.5, 0.05, true, model);
    improvedDifference = abs(tb100 - tb20);
    verifyLessThan(testCase, improvedDifference, 0.2);
    verifyLessThan(testCase, improvedDifference, ...
        abs(base100 - base20) / 20);
    verifyEqual(testCase, detail.radiativeTransferMode, "upperReflection");
end

function reflectivity = local_fresnel_power( ...
        upperPermittivity, lowerPermittivity, sinTheta0, polarization)
    qUpper = sqrt(upperPermittivity - sinTheta0^2);
    qLower = sqrt(lowerPermittivity - sinTheta0^2);
    if polarization == 1
        coefficient = (qUpper - qLower) / (qUpper + qLower);
    else
        coefficient = (lowerPermittivity * qUpper - ...
            upperPermittivity * qLower) / ...
            (lowerPermittivity * qUpper + upperPermittivity * qLower);
    end
    reflectivity = abs(coefficient)^2;
end
