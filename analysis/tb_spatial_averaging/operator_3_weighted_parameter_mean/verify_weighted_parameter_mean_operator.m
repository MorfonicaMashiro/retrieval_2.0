function verify_weighted_parameter_mean_operator()
%VERIFY_WEIGHTED_PARAMETER_MEAN_OPERATOR 验证中心反距离加权参数后算 TB。

    operatorRoot = fileparts(mfilename('fullpath'));
    matlabRoot = fileparts(fileparts(fileparts(operatorRoot)));
    addpath(matlabRoot);

    grid = grid_resampling.ease_grid("EASE");
    coarseRow = 200;
    coarseCol = 200;
    left = grid.xMin + (coarseCol - 1) * grid.cellSize;
    top = grid.yMax - (coarseRow - 1) * grid.cellSize;
    centerX = left + 0.5 * grid.cellSize;
    centerY = top - 0.5 * grid.cellSize;
    x = [centerX + 1000; centerX + 2000];
    y = repmat(centerY, 2, 1);
    [latitude, longitude] = projinv(grid.projection, x, y);

    tb = nan(grid.nRows, grid.nCols);
    tb(coarseRow, coarseCol) = 250;
    smosEase = struct('tb', tb, 'grid', grid);
    thickness = [1; 3];
    snowDepth = [0.1; 0.3];
    oib = make_oib(latitude, longitude, thickness, snowDepth);
    options = struct('computeSampleMeanTB', false, ...
        'parameterWeighting', "inverse_distance_to_center", ...
        'distancePower', 1, 'centerStride', 3, ...
        'centerPhase', [2 2]);
    [pairs, stats] = collocate_day(oib, smosEase, ...
        datetime(2012, 3, 22), struct(), @fake_forward, options);

    expectedArithmeticHi = 2;
    expectedArithmeticHs = 0.2;
    expectedWeightedHi = 5 / 3;
    expectedWeightedHs = 1 / 6;
    assert(height(pairs) == 1 && pairs.nOIB == 2);
    assert(abs(pairs.hi - expectedArithmeticHi) < 1e-12);
    assert(abs(pairs.hs - expectedArithmeticHs) < 1e-12);
    assert(abs(pairs.hiWeighted - expectedWeightedHi) < 1e-12);
    assert(abs(pairs.hsWeighted - expectedWeightedHs) < 1e-12);
    assert(pairs.nOIBWeightedSupport == 2);
    assert(isnan(pairs.nOIBSubgrid));
    assert(abs(pairs.FBiceWeighted - ...
        (pairs.FBsnowWeighted - pairs.hsWeighted)) < 1e-12);
    assert(abs(pairs.TBOIBMeanParameters - ...
        fake_forward(expectedArithmeticHi, expectedArithmeticHs, ...
        false, struct())) < 1e-12);
    assert(abs(pairs.TBOIBWeightedParameters - ...
        fake_forward(expectedWeightedHi, expectedWeightedHs, ...
        false, struct())) < 1e-12);
    assert(~ismember('TBOIBMeanSamples', pairs.Properties.VariableNames), ...
        '算子 3 不应触发逐 OIB 点 TB 计算。');
    assert(stats.forwardPoints == 0 && stats.forwardPairs == 2);
    assert(stats.parameterWeighting == "inverse_distance_to_center");
    assert(stats.distancePower == 1);
    assert(stats.centerStride == 3 && isequal(stats.centerPhase, [2 2]));

    identical = make_oib(latitude, longitude, ...
        repmat(2, 2, 1), repmat(0.2, 2, 1));
    identicalPairs = collocate_day(identical, smosEase, ...
        datetime(2012, 3, 22), struct(), @fake_forward, options);
    assert(abs(identicalPairs.TBWeightedMinusArithmetic) < 1e-12, ...
        'hi/hs 完全相同时，加权和算术参数算子应严格一致。');

    fprintf('verify_weighted_parameter_mean_operator: PASS\n');
end

function oib = make_oib(latitude, longitude, thickness, snowDepth)
    n = numel(latitude);
    meanFreeboard = snowDepth + 0.1;
    tbl = table(latitude, longitude, meanFreeboard, snowDepth, thickness, ...
        zeros(n, 1), zeros(n, 1), zeros(n, 1), zeros(n, 1), ...
        'VariableNames', {'lat','lon','mean_fb','snow_depth','thickness', ...
        'pcnt_ow','pcnt_thin_ice','pcnt_grey_ice','my_ice_flag'});
    oib = struct('source', "synthetic_oib", 'table', tbl);
end

function tb = fake_forward(hi, hs, isMYI, ~)
    tb = hi .^ 2 + 2 .* hs + 10 .* double(isMYI);
end
