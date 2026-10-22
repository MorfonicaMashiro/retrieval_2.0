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
