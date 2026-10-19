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
