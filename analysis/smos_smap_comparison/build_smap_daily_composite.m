function daily = build_smap_daily_composite(filePaths)
%BUILD_SMAP_DAILY_COMPOSITE 将 SPL1CTB 半轨合成为北极 36 km 日亮温。
%   H/V 分别对全部有效前视、后视、升轨和降轨观测作算术平均。

    filePaths = string(filePaths(:));
    if isempty(filePaths)
        error('build_smap_daily_composite:NoFiles', ...
            'At least one SMAP half-orbit file is required.');
    end
    n = 500;
    cellSize = 36000;
    edge = 9000000;
    nCells = n * n;
    sumH = zeros(nCells,1);
    sumV = zeros(nCells,1);
    count = zeros(nCells,1);
    projection = projcrs(6931);
    nObservationsRead = 0;

    for fileIndex = 1:numel(filePaths)
        obs = read_smap_l1c_half_orbit(filePaths(fileIndex));
        nObservationsRead = nObservationsRead + height(obs);
        [x, y] = projfwd(projection, obs.latitude, obs.longitude);
        row = floor((edge - y) ./ cellSize) + 1;
        col = floor((x + edge) ./ cellSize) + 1;
        inside = row >= 1 & row <= n & col >= 1 & col <= n & ...
            isfinite(row) & isfinite(col);
        linear = sub2ind([n n], row(inside), col(inside));
        sumH = sumH + accumarray(linear, obs.tbH(inside), ...
            [nCells 1], @sum, 0);
        sumV = sumV + accumarray(linear, obs.tbV(inside), ...
            [nCells 1], @sum, 0);
        count = count + accumarray(linear, 1, [nCells 1], @sum, 0);
    end

    tbH = reshape(sumH ./ count, [n n]);
    tbV = reshape(sumV ./ count, [n n]);
    observationCount = reshape(count, [n n]);
    tbH(observationCount == 0) = NaN;
    tbV(observationCount == 0) = NaN;
    tbMeanHV = (tbH + tbV) ./ 2;
    xCenters = -edge + ((1:n) - 0.5) * cellSize;
    yCenters = edge - ((1:n) - 0.5) * cellSize;
    daily = struct('tbH', tbH, 'tbV', tbV, 'tbMeanHV', tbMeanHV, ...
        'observationCount', observationCount, 'nRows', n, 'nCols', n, ...
        'cellSize', cellSize, 'edge', edge, 'epsg', 6931, ...
        'projection', projection, 'xCenters', xCenters, ...
        'yCenters', yCenters, 'files', filePaths, ...
        'nObservationsRead', nObservationsRead);
end
