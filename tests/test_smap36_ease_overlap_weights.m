function tests = test_smap36_ease_overlap_weights
%TEST_SMAP36_EASE_OVERLAP_WEIGHTS SMAP36/原版EASE面积交叠回归测试。
    tests = functiontests(localfunctions);
end

function testAreaWeightsCloseForRepresentativeCells(testCase)
    options = struct('edgeSegments', 16, ...
        'areaClosureTolerance', 5e-5);
    [weights, diagnostics] = compute_smap36_ease_overlap_weights( ...
        [200; 250], [300; 250], options);

    verifyEqual(testCase, height(diagnostics), 2);
    verifyLessThan(testCase, ...
        max(abs(diagnostics.areaClosureError)), 5e-5);
    verifyGreaterThan(testCase, min(weights.overlapAreaM2), 0);
    verifyGreaterThan(testCase, min(diagnostics.nEaseCells), 9);
    verifyLessThanOrEqual(testCase, ...
        max(weights.smapAreaFraction), 1);

    for index = 1:height(diagnostics)
        selected = weights.smapRow == diagnostics.smapRow(index) & ...
            weights.smapCol == diagnostics.smapCol(index);
        verifyEqual(testCase, sum(weights.smapAreaFraction(selected)), ...
            1, 'AbsTol', 5e-5);
    end
end

function testDuplicateTargetCellIsComputedOnce(testCase)
    [~, diagnostics] = compute_smap36_ease_overlap_weights( ...
        [250; 250], [250; 250], struct('edgeSegments', 8));
    verifyEqual(testCase, height(diagnostics), 1);
end
