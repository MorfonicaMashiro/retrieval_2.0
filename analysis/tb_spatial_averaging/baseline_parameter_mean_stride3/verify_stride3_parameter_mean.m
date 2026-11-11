function verify_stride3_parameter_mean()
%VERIFY_STRIDE3_PARAMETER_MEAN 验证步长 3 只保留一个非重叠窗口相位。

    experimentRoot = fileparts(mfilename('fullpath'));
    matlabRoot = fileparts(fileparts(fileparts(experimentRoot)));
    addpath(matlabRoot);

    grid = grid_resampling.ease_grid("EASE");
    latitude = repmat(80, 3, 1);
    longitude = zeros(3, 1);
    [pointRow, pointCol] = grid_resampling.PSP2EASE( ...
        latitude(1), longitude(1), grid.version);
    tb = repmat(250, grid.nRows, grid.nCols);
    smosEase = struct('tb', tb, 'grid', grid);
    oib = make_oib(latitude, longitude, [1; 2; 3], [0.1; 0.2; 0.3]);

    allCenterOptions = struct('centerStride', 1, 'centerPhase', [1 1]);
    allCenters = collocate_day( ...
        oib, smosEase, datetime(2012, 3, 22), [], [], allCenterOptions);
    options = struct('centerStride', 3, 'centerPhase', [2 2]);
    [strideCenters, stats] = collocate_day( ...
        oib, smosEase, datetime(2012, 3, 22), [], [], options);

    expected = allCenters( ...
        mod(allCenters.row - 2, 3) == 0 & ...
        mod(allCenters.col - 2, 3) == 0, :);
    assert(height(allCenters) == 9, ...
        '单个内部 OIB 网格应贡献给周围 9 个中心格。');
    assert(height(strideCenters) == 1 && isequaln(strideCenters, expected), ...
        '步长 3 的结果应严格等于全中心格结果的 [2,2] 相位子集。');
    assert(abs(strideCenters.hi - 2) < 1e-12 && ...
        abs(strideCenters.hs - 0.2) < 1e-12, ...
        '步长改变不应改变窗口内原始点算术均值。');
    assert(stats.centerStride == 3 && isequal(stats.centerPhase, [2 2]));
    assert(abs(strideCenters.row - pointRow) <= 1 && ...
        abs(strideCenters.col - pointCol) <= 1);
    fprintf('verify_stride3_parameter_mean: PASS\n');
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
