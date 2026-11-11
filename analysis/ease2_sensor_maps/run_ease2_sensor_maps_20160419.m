function outputs = run_ease2_sensor_maps_20160419(includeSmap9)
%RUN_EASE2_SENSOR_MAPS_20160419 Draw SMOS/SMAP EASE2-N maps.
%   By default, draws SMOS and SMAP 36 km only. Set INCLUDESMAP9=true to
%   add SMAP 9 km when those files are available. All panels use the same
%   2016-04-19 calendar-day files, EASE2-N12.5 km analysis grid,
%   latitude >= 50 deg mask, and shared colour scale. SMAP shows (H+V)/2.

    if nargin < 1 || isempty(includeSmap9), includeSmap9 = false; end
    includeSmap9 = logical(includeSmap9);
    if ~isscalar(includeSmap9)
        error('run_ease2_sensor_maps_20160419:IncludeSMAP9', ...
            'includeSmap9 must be a logical scalar.');
    end

    analysisRoot = fileparts(mfilename('fullpath'));
    matlabRoot = fileparts(fileparts(analysisRoot));
    workspaceRoot = fileparts(matlabRoot);
    dataRoot = fullfile(workspaceRoot, 'data');
    figureRoot = fullfile(analysisRoot, 'figures');
    resultRoot = fullfile(analysisRoot, 'results');
    addpath(matlabRoot, fullfile(matlabRoot, 'readdata'), ...
        fullfile(matlabRoot, 'readdata', 'readsmap'), ...
        fullfile(matlabRoot, 'analysis', 'smos_smap_comparison'), ...
        analysisRoot);
    if ~isfolder(figureRoot), mkdir(figureRoot); end
    if ~isfolder(resultRoot), mkdir(resultRoot); end

    dayToken = '20160419';
    smosFile = fullfile(dataRoot, 'smos', '2016', ...
        ['SMOS_TB_north_' dayToken '.nc']);
    smap9Directory = fullfile(dataRoot, 'smap', 'SPL1CTB_E_004', ...
        '2016', '04', '19');
    smap36Directory = fullfile(dataRoot, 'smap', 'SPL1CTB_006', ...
        '2016', '04', '19');
    require_file(smosFile, 'SMOS');
    smap36Files = list_calendar_day_files(smap36Directory, ...
        ['SMAP_L1C_TB_*_' dayToken 'T*_R*_*.h5']);
    if numel(smap36Files) ~= 29
        error('run_ease2_sensor_maps_20160419:SMAP36Count', ...
            'Expected 29 SMAP 36 km files for %s; found %d.', ...
            dayToken, numel(smap36Files));
    end

    smos = grid_resampling.regrid_smos_to_ease( ...
        read_smos(smosFile), "EASE2");
    smap36 = build_smap_daily_ease2_12p5(smap36Files, 36);
    grid = smos.grid;
    latitudeMask = build_latitude_mask(grid, 50);
    smosTB = mask_field(smos.tb, latitudeMask);
    smap36TB = mask_field(smap36.tbMeanHV, latitudeMask);
    maps = {smosTB, smap36TB};
    titles = {'SMOS L3B daily TB', ...
        'SMAP standard 36 km daily (H+V)/2'};
    subtitles = {'nearest-neighbour remap from EPSG:3413', ...
        'native EASE2-N 36 km; sampled at 12.5 km centres'};
    stems = {'smos', 'smap36km'};
    smap9Files = strings(0, 1);
    if includeSmap9
        smap9Files = list_calendar_day_files(smap9Directory, ...
            ['SMAP_L1C_TB_E_*_' dayToken 'T*_R*_*.h5']);
        if numel(smap9Files) ~= 29
            error('run_ease2_sensor_maps_20160419:SMAP9Count', ...
                'Expected 29 SMAP 9 km files for %s; found %d.', ...
                dayToken, numel(smap9Files));
        end
        smap9 = build_smap_daily_ease2_12p5(smap9Files, 9);
        smap9TB = mask_field(smap9.tbMeanHV, latitudeMask);
        maps = {smosTB, smap9TB, smap36TB};
        titles = {'SMOS L3B daily TB', ...
            'SMAP enhanced 9 km daily (H+V)/2', ...
            'SMAP standard 36 km daily (H+V)/2'};
        subtitles = {'nearest-neighbour remap from EPSG:3413', ...
            'native EASE2-N 9 km; sampled at 12.5 km centres', ...
            'native EASE2-N 36 km; sampled at 12.5 km centres'};
        stems = {'smos', 'smap9km', 'smap36km'};
        outputStem = 'smos_smap9_smap36_ease2_20160419';
    else
        outputStem = 'smos_smap36_ease2_20160419';
    end

    values = cellfun(@(map) map(isfinite(map)), maps, ...
        'UniformOutput', false);
    values = vertcat(values{:});
    commonLimits = [floor(prctile(values, 0.5) / 5) * 5, ...
        ceil(prctile(values, 99.5) / 5) * 5];
    commonLimits(1) = max(50, commonLimits(1));
    commonLimits(2) = min(310, commonLimits(2));
    nMaps = numel(maps);

    combinedPath = fullfile(figureRoot, [outputStem '.png']);
    combinedPdf = "";
    pdfCreated = false;
    fig = figure('Visible', 'off', 'Color', 'w', ...
        'Position', [40 80 620 * nMaps 620]);
    layout = tiledlayout(fig, 1, nMaps, ...
        'TileSpacing', 'compact', 'Padding', 'compact');
    title(layout, ...
        '2016-04-19 brightness temperature on EASE-Grid 2.0 North', ...
        'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 15);
    axesList = gobjects(1, nMaps);
    for index = 1:nMaps
        axesList(index) = nexttile(layout);
        draw_map(axesList(index), maps{index}, grid, ...
            titles{index}, subtitles{index}, commonLimits);
    end
    colormap(fig, parula(256));
    colorbarHandle = colorbar(axesList(end), 'eastoutside');
    colorbarHandle.Label.String = 'Brightness temperature (K)';
    colorbarHandle.Label.FontName = 'Arial';
    colorbarHandle.FontName = 'Arial';
    exportgraphics(fig, combinedPath, 'Resolution', 300);
    close(fig);

    individualPaths = strings(nMaps, 1);
    for index = 1:nMaps
        individualPaths(index) = fullfile(figureRoot, ...
            sprintf('%s_ease2_20160419.png', stems{index}));
        fig = figure('Visible', 'off', 'Color', 'w', ...
            'Position', [80 80 760 670]);
        ax = axes(fig); %#ok<LAXES>
        draw_map(ax, maps{index}, grid, titles{index}, ...
            subtitles{index}, commonLimits);
        colormap(fig, parula(256));
        colorbarHandle = colorbar(ax, 'eastoutside');
        colorbarHandle.Label.String = 'Brightness temperature (K)';
        colorbarHandle.FontName = 'Arial';
        colorbarHandle.Label.FontName = 'Arial';
        exportgraphics(fig, individualPaths(index), 'Resolution', 300);
        close(fig);
    end

    statistics = build_statistics(maps, titles, commonLimits);
    if includeSmap9
        statistics.nFiles = [1; numel(smap9Files); numel(smap36Files)];
    else
        statistics.nFiles = [1; numel(smap36Files)];
    end
    writetable(statistics, fullfile(resultRoot, ...
        'ease2_maps_20160419_statistics.csv'));
    save(fullfile(resultRoot, [outputStem '.mat']), ...
        'maps', 'titles', 'grid', 'commonLimits', 'statistics', '-v7.3');

    outputs = struct('combinedPng', string(combinedPath), ...
        'combinedPdf', combinedPdf, ...
        'pdfCreated', pdfCreated, ...
        'individualPng', individualPaths, ...
        'statistics', statistics, 'commonLimitsK', commonLimits, ...
        'smap9FileCount', numel(smap9Files), ...
        'smap36FileCount', numel(smap36Files));
    fprintf('Created EASE2 sensor maps:\n%s\n', combinedPath);
end

function files = list_calendar_day_files(directory, pattern)
    if ~isfolder(directory)
        error('run_ease2_sensor_maps_20160419:MissingDirectory', ...
            'Directory does not exist: %s', directory);
    end
    listing = dir(fullfile(directory, pattern));
    [~, order] = sort({listing.name});
    listing = listing(order);
    files = string(fullfile({listing.folder}, {listing.name})).';
end

function require_file(filePath, label)
    if ~isfile(filePath)
        error('run_ease2_sensor_maps_20160419:MissingFile', ...
            '%s file not found: %s', label, filePath);
    end
end

function mask = build_latitude_mask(grid, minimumLatitude)
    [xLimit, yLimit] = projfwd(grid.projection, minimumLatitude, 0);
    radius = hypot(xLimit, yLimit);
    mask = grid.yCenters(:).^2 + grid.xCenters(:)'.^2 <= radius^2;
end

function field = mask_field(field, mask)
    field = double(field);
    field(~mask) = NaN;
end

function draw_map(ax, field, grid, mainTitle, subTitle, colorLimits)
    xKm = grid.xCenters / 1000;
    yKm = fliplr(grid.yCenters) / 1000;
    imageHandle = imagesc(ax, xKm, yKm, flipud(field));
    imageHandle.AlphaData = isfinite(flipud(field));
    ax.Color = [0.92 0.94 0.96];
    ax.YDir = 'normal';
    axis(ax, 'image');
    xlim(ax, [-4700 4700]);
    ylim(ax, [-4700 4700]);
    clim(ax, colorLimits);
    hold(ax, 'on');
    draw_graticule(ax, grid.projection);
    draw_coastline(ax, grid.projection);
    hold(ax, 'off');
    title(ax, {mainTitle; subTitle}, 'Interpreter', 'none', ...
        'FontName', 'Arial', 'FontWeight', 'normal', 'FontSize', 10.5);
    xlabel(ax, 'EASE2 x (km)');
    ylabel(ax, 'EASE2 y (km)');
    ax.FontName = 'Arial';
    ax.FontSize = 9;
    ax.LineWidth = 0.8;
    ax.TickDir = 'out';
    box(ax, 'on');
end

function draw_graticule(ax, projection)
    longitude = linspace(-180, 180, 721);
    for latitude = [50 60 70 80]
        [x, y] = projfwd(projection, ...
            latitude * ones(size(longitude)), longitude);
        plot(ax, x / 1000, y / 1000, '-', ...
            'Color', [0.88 0.88 0.88], 'LineWidth', 0.55);
    end
    latitude = linspace(50, 89.9, 240);
    for longitudeValue = -180:45:135
        [x, y] = projfwd(projection, latitude, ...
            longitudeValue * ones(size(latitude)));
        plot(ax, x / 1000, y / 1000, '-', ...
            'Color', [0.88 0.88 0.88], 'LineWidth', 0.5);
    end
end

function draw_coastline(ax, projection)
    coast = load('coastlines');
    valid = coast.coastlat >= 49 | isnan(coast.coastlat);
    latitude = coast.coastlat;
    longitude = coast.coastlon;
    latitude(~valid) = NaN;
    longitude(~valid) = NaN;
    [x, y] = projfwd(projection, latitude, longitude);
    plot(ax, x / 1000, y / 1000, '-', ...
        'Color', [0.12 0.12 0.12], 'LineWidth', 0.55);
end

function statistics = build_statistics(maps, names, colorLimits)
    nMaps = numel(maps);
    n = zeros(nMaps, 1);
    minimum = nan(nMaps, 1);
    medianValue = nan(nMaps, 1);
    maximum = nan(nMaps, 1);
    p01 = nan(nMaps, 1);
    p99 = nan(nMaps, 1);
    for index = 1:nMaps
        values = maps{index}(isfinite(maps{index}));
        n(index) = numel(values);
        minimum(index) = min(values);
        medianValue(index) = median(values);
        maximum(index) = max(values);
        p01(index) = prctile(values, 1);
        p99(index) = prctile(values, 99);
    end
    product = string(names(:));
    colorMinK = repmat(colorLimits(1), nMaps, 1);
    colorMaxK = repmat(colorLimits(2), nMaps, 1);
    statistics = table(product, n, minimum, p01, medianValue, p99, ...
        maximum, colorMinK, colorMaxK);
end
