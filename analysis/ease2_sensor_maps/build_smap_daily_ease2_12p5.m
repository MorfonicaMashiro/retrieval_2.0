function daily = build_smap_daily_ease2_12p5(filePaths, nativeSpacingKm)
%BUILD_SMAP_DAILY_EASE2_12P5 SMAP half-orbits to a daily EASE2-N12.5 raster.
%   DAILY = BUILD_SMAP_DAILY_EASE2_12P5(FILEPATHS, NATIVESPACINGKM)
%   averages valid fore/aft H and V observations on the native EASE2-N grid,
%   then samples that daily native raster at EASE2-N12.5 km cell centres.
%
%   The 12.5 km raster is an analysis grid only. Nearest-neighbour sampling
%   does not increase the native 9 km or 36 km physical resolution.

    filePaths = string(filePaths(:));
    nativeSpacingKm = double(nativeSpacingKm);
    if isempty(filePaths)
        error('build_smap_daily_ease2_12p5:NoFiles', ...
            'At least one SMAP half-orbit file is required.');
    end
    if ~isscalar(nativeSpacingKm) || ...
            ~ismember(nativeSpacingKm, [9, 36])
        error('build_smap_daily_ease2_12p5:Spacing', ...
            'nativeSpacingKm must be 9 or 36.');
    end
    missing = filePaths(~isfile(filePaths));
    if ~isempty(missing)
        error('build_smap_daily_ease2_12p5:MissingFile', ...
            'SMAP file does not exist: %s', missing(1));
    end

    nativeCellSize = nativeSpacingKm * 1000;
    edge = 9000000;
    nativeSize = round(2 * edge / nativeCellSize);
    if nativeSize * nativeCellSize ~= 2 * edge
        error('build_smap_daily_ease2_12p5:Geometry', ...
            'Native spacing does not divide the EASE2-N extent.');
    end
    nNativeCells = nativeSize ^ 2;
    sumH = zeros(nNativeCells, 1);
    sumV = zeros(nNativeCells, 1);
    count = zeros(nNativeCells, 1, 'uint32');
    projection = projcrs(6931);
    observationsRead = 0;
    observationsUsed = 0;

    for fileIndex = 1:numel(filePaths)
        observations = read_smap_l1c_half_orbit(filePaths(fileIndex));
        observationsRead = observationsRead + height(observations);
        valid = observations.tbH <= 310 & observations.tbV <= 310;
        observations = observations(valid, :);
        [x, y] = projfwd(projection, ...
            observations.latitude, observations.longitude);
        row = floor((edge - y) ./ nativeCellSize) + 1;
        col = floor((x + edge) ./ nativeCellSize) + 1;
        inside = row >= 1 & row <= nativeSize & ...
            col >= 1 & col <= nativeSize & ...
            isfinite(row) & isfinite(col);
        if ~any(inside)
            continue
        end
        linear = sub2ind([nativeSize nativeSize], ...
            row(inside), col(inside));
        sumH = sumH + accumarray(linear, observations.tbH(inside), ...
            [nNativeCells 1], @sum, 0);
        sumV = sumV + accumarray(linear, observations.tbV(inside), ...
            [nNativeCells 1], @sum, 0);
        count = count + uint32(accumarray(linear, 1, ...
            [nNativeCells 1], @sum, 0));
        observationsUsed = observationsUsed + nnz(inside);
    end

    nativeCount = reshape(count, [nativeSize nativeSize]);
    nativeH = reshape(sumH, [nativeSize nativeSize]) ./ ...
        double(nativeCount);
    nativeV = reshape(sumV, [nativeSize nativeSize]) ./ ...
        double(nativeCount);
    nativeH(nativeCount == 0) = NaN;
    nativeV(nativeCount == 0) = NaN;
    nativeMeanHV = (nativeH + nativeV) ./ 2;

    grid = grid_resampling.ease_grid("EASE2");
    nativeX0 = -edge + nativeCellSize / 2;
    nativeY0 = edge - nativeCellSize / 2;
    sourceCol = round((grid.xCenters - nativeX0) ./ nativeCellSize) + 1;
    sourceRow = round((nativeY0 - grid.yCenters) ./ nativeCellSize) + 1;
    sourceCol = min(max(sourceCol, 1), nativeSize);
    sourceRow = min(max(sourceRow, 1), nativeSize);

    daily = struct();
    daily.tbH = nativeH(sourceRow, sourceCol);
    daily.tbV = nativeV(sourceRow, sourceCol);
    daily.tbMeanHV = nativeMeanHV(sourceRow, sourceCol);
    daily.observationCount = nativeCount(sourceRow, sourceCol);
    daily.grid = grid;
    daily.nativeSpacingKm = nativeSpacingKm;
    daily.nativeSize = nativeSize;
    daily.nativeTBMeanHV = nativeMeanHV;
    daily.nativeObservationCount = nativeCount;
    daily.files = filePaths;
    daily.nFiles = numel(filePaths);
    daily.nObservationsRead = observationsRead;
    daily.nObservationsUsed = observationsUsed;
    daily.method = "native daily mean, then nearest native cell at target centre";
    daily.resolutionNote = ...
        "12.5 km is target spacing; native physical resolution is unchanged.";
end
