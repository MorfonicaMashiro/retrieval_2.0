function tests = test_regrid_to_ease2_12p5
%TEST_REGRID_TO_EASE2_12P5 统一 EASE2-N12.5km 转换入口回归测试。
    tests = functiontests(localfunctions);
end

function testOfficialGridDefinition(testCase)
    grid = grid_resampling.ease_grid("EASE2");
    verifyEqual(testCase, grid.epsg, 6931);
    verifyEqual(testCase, [grid.nRows, grid.nCols], [1440, 1440]);
    verifyEqual(testCase, grid.cellSize, 12500);
    verifyEqual(testCase, grid.xMin, -9000000);
    verifyEqual(testCase, grid.yMax, 9000000);
    verifyEqual(testCase, grid.xCenters([1 end]), [-8993750, 8993750]);
    verifyEqual(testCase, grid.yCenters([1 end]), [8993750, -8993750]);
end

function testOibRecordsAreLocatedWithoutDataLoss(testCase)
    grid = grid_resampling.ease_grid("EASE2");
    targetRow = 100;
    targetCol = 200;
    [latitude, longitude] = projinv(grid.projection, ...
        grid.xCenters(targetCol), grid.yCenters(targetRow));
    lat = [latitude; latitude; NaN];
    lon = [longitude; longitude; 0];
    thickness = [1.1; 1.2; 1.3];
    oib = struct('source', "fixture.txt", ...
        'table', table(lat, lon, thickness));

    result = grid_resampling.regrid_to_ease2_12p5(oib);

    verifyEqual(testCase, height(result.records), 3);
    verifyEqual(testCase, height(result.validRecords), 2);
    verifyEqual(testCase, result.validRecords.ease2Row, ...
        repmat(targetRow, 2, 1));
    verifyEqual(testCase, result.validRecords.ease2Col, ...
        repmat(targetCol, 2, 1));
    verifyEqual(testCase, result.count(targetRow, targetCol), uint32(2));
    verifyEqual(testCase, result.records.thickness, thickness);
    verifyFalse(testCase, result.physicalResolutionChanged);
end

function testSmapSelectionPreservesNativeResolution(testCase)
    grid = grid_resampling.ease_grid("EASE2");
    targetRow = 700;
    targetCol = 710;
    [latitude, longitude] = projinv(grid.projection, ...
        grid.xCenters(targetCol), grid.yCenters(targetRow));
    observations = table( ...
        [latitude; latitude], [longitude; longitude], ...
        [true; false], [250; 251], ...
        'VariableNames', {'latitude','longitude','recommendedMask','tbH'});
    data = struct('source', "fixture.h5", ...
        'product', "SPL1CTB.006", 'resolutionKm', 36, ...
        'observations', observations, ...
        'recommended', observations(observations.recommendedMask, :));
    options = struct('smapSelection', "recommended");

    result = grid_resampling.regrid_to_ease2_12p5( ...
        data, "smap", options);

    verifyEqual(testCase, height(result.records), 1);
    verifyEqual(testCase, result.sourceResolutionKm, 36);
    verifyEqual(testCase, result.targetSpacingKm, 12.5);
    verifyEqual(testCase, result.count(targetRow, targetCol), uint32(1));
    verifyEqual(testCase, result.method, ...
        "point_to_cell_without_time_aggregation");
end

function testUnknownInputIsRejected(testCase)
    verifyError(testCase, ...
        @() grid_resampling.regrid_to_ease2_12p5(struct()), ...
        'regrid_to_ease2_12p5:AutoDetect');
end
