function [weights, diagnostics] = compute_smap36_ease_overlap_weights( ...
        smapRow, smapCol, options)
%COMPUTE_SMAP36_EASE_OVERLAP_WEIGHTS 计算 SMAP36 与原版 EASE 的面积交叠。
%   [WEIGHTS,DIAGNOSTICS] = COMPUTE_SMAP36_EASE_OVERLAP_WEIGHTS(ROW,COL)
%   把每个 36 km EASE2-N SMAP 网格作为目标支撑域，将原版北极
%   EASE 12.53376 km 像元边界转换到 EASE2 等面积坐标后作多边形相交。
%   WEIGHTS 中的 smapAreaFraction 是交叠面积占 36 km 目标格面积的比例。
%
%   为表示投影转换后的曲线边界，每条原版 EASE 像元边默认细分为
%   32 段。该方法是面积多边形相交，不是格心归属或固定 3x3 近似。

    if nargin < 3 || isempty(options), options = struct(); end
    options = validate_options(options);
    smapRow = double(smapRow(:));
    smapCol = double(smapCol(:));
    if numel(smapRow) ~= numel(smapCol) || isempty(smapRow) || ...
            any(~isfinite(smapRow)) || any(~isfinite(smapCol)) || ...
            any(smapRow ~= fix(smapRow)) || any(smapCol ~= fix(smapCol))
        error('compute_smap36_ease_overlap_weights:Cells', ...
            'smapRow and smapCol must be nonempty equal-length integer vectors.');
    end

    smapGrid = struct('nRows', 500, 'nCols', 500, ...
        'cellSize', 36000, 'xMin', -9000000, 'yMax', 9000000, ...
        'projection', projcrs(6931));
    if any(smapRow < 1 | smapRow > smapGrid.nRows | ...
            smapCol < 1 | smapCol > smapGrid.nCols)
        error('compute_smap36_ease_overlap_weights:Cells', ...
            'SMAP 36 km row/column must be within 1..500.');
    end
    [cells, ~] = unique([smapRow smapCol], 'rows', 'stable');
    smapRow = cells(:,1);
    smapCol = cells(:,2);

    sourceGrid = grid_resampling.ease_grid("EASE");
    targetAreaM2 = smapGrid.cellSize ^ 2;
    parts = cell(size(smapRow));
    diagnosticParts = cell(size(smapRow));

    for cellIndex = 1:numel(smapRow)
        row = smapRow(cellIndex);
        col = smapCol(cellIndex);
        xLeft = smapGrid.xMin + (col - 1) * smapGrid.cellSize;
        xRight = xLeft + smapGrid.cellSize;
        yTop = smapGrid.yMax - (row - 1) * smapGrid.cellSize;
        yBottom = yTop - smapGrid.cellSize;
        targetPolygon = polyshape( ...
            [xLeft xRight xRight xLeft], ...
            [yBottom yBottom yTop yTop], 'Simplify', false);

        [targetBoundaryX, targetBoundaryY] = rectangle_boundary( ...
            xLeft, xRight, yBottom, yTop, options.edgeSegments);
        [boundaryLatitude, boundaryLongitude] = projinv( ...
            smapGrid.projection, targetBoundaryX, targetBoundaryY);
        [boundarySourceX, boundarySourceY] = projfwd( ...
            sourceGrid.projection, boundaryLatitude, boundaryLongitude);
        candidateCols = candidate_indices( ...
            min(boundarySourceX), max(boundarySourceX), ...
            sourceGrid.xMin, sourceGrid.cellSize, sourceGrid.nCols);
        candidateRows = candidate_rows( ...
            min(boundarySourceY), max(boundarySourceY), ...
            sourceGrid.yMax, sourceGrid.cellSize, sourceGrid.nRows);

        overlapRows = zeros(numel(candidateRows) * numel(candidateCols), 1);
        overlapCols = zeros(size(overlapRows));
        overlapAreas = zeros(size(overlapRows));
        kept = 0;
        for rowIndex = 1:numel(candidateRows)
            sourceRow = candidateRows(rowIndex);
            sourceYTop = sourceGrid.yMax - ...
                (sourceRow - 1) * sourceGrid.cellSize;
            sourceYBottom = sourceYTop - sourceGrid.cellSize;
            for colIndex = 1:numel(candidateCols)
                sourceCol = candidateCols(colIndex);
                sourceXLeft = sourceGrid.xMin + ...
                    (sourceCol - 1) * sourceGrid.cellSize;
                sourceXRight = sourceXLeft + sourceGrid.cellSize;
                [sourceX, sourceY] = rectangle_boundary( ...
                    sourceXLeft, sourceXRight, sourceYBottom, sourceYTop, ...
                    options.edgeSegments);
                [latitude, longitude] = projinv( ...
                    sourceGrid.projection, sourceX, sourceY);
                [targetX, targetY] = projfwd( ...
                    smapGrid.projection, latitude, longitude);
                sourcePolygon = polyshape( ...
                    targetX, targetY, 'Simplify', false);
                overlapPolygon = intersect(targetPolygon, sourcePolygon);
                overlapAreaM2 = area(overlapPolygon);
                if overlapAreaM2 > options.minimumAreaM2
                    kept = kept + 1;
                    overlapRows(kept) = sourceRow;
                    overlapCols(kept) = sourceCol;
                    overlapAreas(kept) = overlapAreaM2;
                end
            end
        end
        overlapRows = overlapRows(1:kept);
        overlapCols = overlapCols(1:kept);
        overlapAreas = overlapAreas(1:kept);
        areaFraction = overlapAreas ./ targetAreaM2;
        parts{cellIndex} = table( ...
            repmat(row, kept, 1), repmat(col, kept, 1), ...
            overlapRows, overlapCols, overlapAreas, ...
            overlapAreas ./ 1e6, areaFraction, ...
            'VariableNames', {'smapRow','smapCol','easeRow','easeCol', ...
            'overlapAreaM2','overlapAreaKm2','smapAreaFraction'});
        diagnosticParts{cellIndex} = table(row, col, kept, ...
            sum(overlapAreas), sum(areaFraction), ...
            'VariableNames', {'smapRow','smapCol','nEaseCells', ...
            'coveredAreaM2','geometryCoverageFraction'});
    end

    weights = vertcat(parts{:});
    diagnostics = vertcat(diagnosticParts{:});
    diagnostics.areaClosureError = ...
        diagnostics.geometryCoverageFraction - 1;
    if any(abs(diagnostics.areaClosureError) > options.areaClosureTolerance)
        bad = diagnostics(abs(diagnostics.areaClosureError) > ...
            options.areaClosureTolerance, :);
        error('compute_smap36_ease_overlap_weights:AreaClosure', ...
            ['Area overlap does not close for %d SMAP cells; maximum ' ...
            'absolute fractional error is %.3g.'], ...
            height(bad), max(abs(bad.areaClosureError)));
    end
end

function options = validate_options(options)
    defaults = struct('edgeSegments', 32, 'minimumAreaM2', 1, ...
        'areaClosureTolerance', 2e-5);
    if ~isstruct(options) || ~isscalar(options)
        error('compute_smap36_ease_overlap_weights:Options', ...
            'options must be a scalar struct.');
    end
    names = fieldnames(options);
    unknown = setdiff(names, fieldnames(defaults));
    if ~isempty(unknown)
        error('compute_smap36_ease_overlap_weights:Options', ...
            'Unknown option: %s.', strjoin(unknown, ', '));
    end
    for index = 1:numel(names)
        defaults.(names{index}) = options.(names{index});
    end
    options = defaults;
    if ~isscalar(options.edgeSegments) || options.edgeSegments < 2 || ...
            options.edgeSegments ~= fix(options.edgeSegments)
        error('compute_smap36_ease_overlap_weights:Options', ...
            'edgeSegments must be an integer >= 2.');
    end
end

function indices = candidate_indices( ...
        minimumX, maximumX, xMinimum, cellSize, nColumns)
    first = floor((minimumX - xMinimum) / cellSize) + 1;
    last = floor((maximumX - xMinimum) / cellSize) + 1;
    indices = max(1, first - 1):min(nColumns, last + 1);
end

function rows = candidate_rows( ...
        minimumY, maximumY, yMaximum, cellSize, nRows)
    first = floor((yMaximum - maximumY) / cellSize) + 1;
    last = floor((yMaximum - minimumY) / cellSize) + 1;
    rows = max(1, first - 1):min(nRows, last + 1);
end

function [x, y] = rectangle_boundary( ...
        xLeft, xRight, yBottom, yTop, segments)
    bottomX = linspace(xLeft, xRight, segments + 1);
    bottomY = repmat(yBottom, 1, segments + 1);
    rightX = repmat(xRight, 1, segments);
    rightY = linspace(yBottom, yTop, segments + 1);
    rightY = rightY(2:end);
    topX = linspace(xRight, xLeft, segments + 1);
    topX = topX(2:end);
    topY = repmat(yTop, 1, segments);
    leftX = repmat(xLeft, 1, segments - 1);
    leftY = linspace(yTop, yBottom, segments + 1);
    leftY = leftY(2:end-1);
    x = [bottomX rightX topX leftX];
    y = [bottomY rightY topY leftY];
end
