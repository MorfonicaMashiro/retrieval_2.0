function verify_sample_tb_mean_operator()
%VERIFY_SAMPLE_TB_MEAN_OPERATOR 验证逐 OIB 点算 TB 后平均的算子。

    operatorRoot = fileparts(mfilename('fullpath'));
    matlabRoot = fileparts(fileparts(fileparts(operatorRoot)));
    addpath(matlabRoot);

    grid = grid_resampling.ease_grid("EASE");
    latitude = repmat(80, 3, 1);
    longitude = zeros(3, 1);
    [row, col] = grid_resampling.PSP2EASE( ...
        latitude(1), longitude(1), grid.version);
    tb = nan(grid.nRows, grid.nCols);
    tb(row, col) = 250;
    smosEase = struct('tb', tb, 'grid', grid);
    model = struct();
    phase = [mod(row - 1, 3) + 1, mod(col - 1, 3) + 1];
    options = struct('centerStride', 3, 'centerPhase', phase);

    varied = make_oib(latitude, longitude, [1; 2; 3], [0.1; 0.2; 0.3]);
    [baseline, ~] = collocate_day( ...
        varied, smosEase, datetime(2012, 3, 22));
    assert(isequal(baseline.Properties.VariableNames, ...
        {'day','flight','row','col','TB','FBsnow','hs','FBice', ...
        'hi','nOIB','isMYI'}), ...
        '三参数调用的历史输出字段发生了变化。');

    [pairs, stats] = collocate_day( ...
        varied, smosEase, datetime(2012, 3, 22), ...
        model, @fake_forward, options);
    expectedMeanParameters = fake_forward(mean([1; 2; 3]), ...
        mean([0.1; 0.2; 0.3]), false, model);
    expectedMeanSamples = mean(fake_forward( ...
        [1; 2; 3], [0.1; 0.2; 0.3], false, model));
    assert(height(pairs) == 1 && pairs.nOIB == 3);
    assert(abs(pairs.TBOIBMeanParameters - expectedMeanParameters) < 1e-12);
    assert(abs(pairs.TBOIBMeanSamples - expectedMeanSamples) < 1e-12);
    assert(abs(pairs.TBSampleMinusParameter - ...
        (expectedMeanSamples - expectedMeanParameters)) < 1e-12);
    assert(stats.forwardPoints == 3 && isfinite(stats.forwardSeconds));

    identical = make_oib(latitude, longitude, ...
        ones(3, 1), repmat(0.1, 3, 1));
    identicalPairs = collocate_day( ...
        identical, smosEase, datetime(2012, 3, 22), ...
        model, @fake_forward, options);
    assert(abs(identicalPairs.TBSampleMinusParameter) < 1e-12, ...
        'hi/hs 完全相同时，两种平均次序应严格一致。');

    fprintf('verify_sample_tb_mean_operator: PASS\n');
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
