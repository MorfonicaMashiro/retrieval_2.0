function [solutions, info] = invert_pointwise_freeboard_pair( ...
        tbObs, fbSamples, branch, isMYI, model, forwardTB, weights, options)
%INVERT_POINTWISE_FREEBOARD_PAIR 用逐点干舷构造足迹平均 TB 并求根。
%   OIB hi/hs 真值不进入本函数。雪面干舷路线采用 2018 论文的平坦雪层
%   扫描：候选雪深超过局部雪面干舷时，将该点雪深限制为其雪面干舷；
%   潜在淹没的加权占比超过阈值后，该候选值无效。

    if nargin < 7 || isempty(weights), weights = ones(size(fbSamples)); end
    if nargin < 8 || isempty(options), options = struct(); end
    defaults = struct('maxSnowDepthM', 1, 'maxIceThicknessM', 10, ...
        'initialStepM', 0.04, 'scanStepM', 0.01, ...
        'tbToleranceK', 0.01, 'snowToleranceM', 1e-5, ...
        'maxIterations', 60, 'maxInundatedFraction', 0.5);
    names = fieldnames(defaults);
    for index = 1:numel(names)
        if ~isfield(options, names{index})
            options.(names{index}) = defaults.(names{index});
        end
    end

    branch = string(branch);
    fbSamples = double(fbSamples(:));
    weights = double(weights(:));
    required = {'rhoWater','rhoIce','rhoSnow'};
    if ~isscalar(tbObs) || ~isfinite(tbObs) || isempty(fbSamples) || ...
            numel(weights) ~= numel(fbSamples) || ...
            any(~isfinite(fbSamples)) || ...
            any(~isfinite(weights) | weights <= 0) || ...
            ~ismember(branch, ["ice","snow"]) || ...
            ~isscalar(isMYI) || ~ismember(double(isMYI), [0 1]) || ...
            ~isstruct(model) || ~all(isfield(model, required)) || ...
            ~isa(forwardTB, 'function_handle')
        error('invert_pointwise_freeboard_pair:Inputs', ...
            'TB、逐点干舷、权重、冰型、密度或前向函数无效。');
    end
    weights = weights / sum(weights);
    rhoW = model.rhoWater; rhoI = model.rhoIce; rhoS = model.rhoSnow;
    if ~(rhoW > rhoI && rhoI > rhoS && rhoS > 0) || ...
            options.maxSnowDepthM <= 0 || options.maxIceThicknessM <= 0 || ...
            options.initialStepM <= 0 || options.scanStepM <= 0 || ...
            options.tbToleranceK <= 0 || options.snowToleranceM <= 0 || ...
            options.maxIterations < 1 || options.maxInundatedFraction < 0 || ...
            options.maxInundatedFraction > 1
        error('invert_pointwise_freeboard_pair:Parameters', ...
            '密度或求解配置不合理。');
    end

    solutions = struct('hi', {}, 'hs', {}, 'candidateHs', {}, ...
        'tbSim', {}, 'tbResidual', {}, 'inundatedFraction', {}, ...
        'nSamples', {});
    info = struct('status', "no_root", 'scanCount', 0, 'rootCount', 0, ...
        'unconvergedBrackets', 0, 'scanTBMinK', NaN, ...
        'scanTBMaxK', NaN);
    if branch == "snow" && any(fbSamples < 0)
        info.status = "negative_snow_freeboard";
        return
    end

    step = options.scanStepM;
    if branch == "ice", step = options.initialStepM; end
    hsGrid = 0:step:options.maxSnowDepthM;
    if hsGrid(end) < options.maxSnowDepthM
        hsGrid(end + 1) = options.maxSnowDepthM;
    end
    residual = nan(size(hsGrid));
    for gridIndex = 1:numel(hsGrid)
        [~, ~, ~, residual(gridIndex)] = evaluate(hsGrid(gridIndex));
    end
    info.scanCount = numel(hsGrid);
    finiteResidual = residual(isfinite(residual));
    if ~isempty(finiteResidual)
        info.scanTBMinK = tbObs + min(finiteResidual);
        info.scanTBMaxK = tbObs + max(finiteResidual);
    end

    hit = isfinite(residual) & abs(residual) <= options.tbToleranceK;
    edges = diff([false, hit, false]);
    first = find(edges == 1);
    last = find(edges == -1) - 1;
    for hitIndex = 1:numel(first)
        add_solution(hsGrid(round((first(hitIndex) + last(hitIndex)) / 2)));
    end
    for gridIndex = 1:numel(hsGrid)-1
        if hit(gridIndex) || hit(gridIndex+1) || ...
                ~isfinite(residual(gridIndex)) || ...
                ~isfinite(residual(gridIndex+1)) || ...
                residual(gridIndex) * residual(gridIndex+1) >= 0
            continue
        end
        lo = hsGrid(gridIndex); hiBracket = hsGrid(gridIndex+1);
        rLo = residual(gridIndex); converged = false;
        for iteration = 1:options.maxIterations
            middle = (lo + hiBracket) / 2;
            [~, ~, ~, rMiddle] = evaluate(middle);
            if ~isfinite(rMiddle), break; end
            if abs(rMiddle) <= options.tbToleranceK
                converged = add_solution(middle);
                break
            end
            if rLo * rMiddle <= 0
                hiBracket = middle;
            else
                lo = middle; rLo = rMiddle;
            end
            if hiBracket - lo <= options.snowToleranceM
                converged = add_solution((lo + hiBracket) / 2);
                break
            end
        end
        if ~converged
            info.unconvergedBrackets = info.unconvergedBrackets + 1;
        end
    end
    if ~isempty(solutions)
        [~, order] = sort([solutions.candidateHs]);
        solutions = solutions(order);
        info.rootCount = numel(solutions);
        info.status = "solved";
        if branch == "ice" && info.rootCount > 1
            info.status = "multiple_radar_roots";
        end
    end

    function [hiMean, hsMean, tbMean, r, inundatedFraction] = evaluate(hsValue)
        hiMean = NaN; hsMean = NaN; tbMean = NaN; r = NaN;
        inundatedFraction = 0;
        if branch == "ice"
            hsLocal = repmat(hsValue, size(fbSamples));
            hiLocal = (rhoW .* fbSamples + rhoS .* hsLocal) ./ (rhoW-rhoI);
        else
            inundated = fbSamples < hsValue;
            inundatedFraction = sum(weights(inundated));
            if inundatedFraction > options.maxInundatedFraction, return; end
            hsLocal = min(hsValue, fbSamples);
            hiLocal = (rhoW .* fbSamples - (rhoW-rhoS) .* hsLocal) ...
                ./ (rhoW-rhoI);
        end
        if any(~isfinite(hiLocal) | ~isfinite(hsLocal) | ...
                hiLocal < 0 | hiLocal > options.maxIceThicknessM | hsLocal < 0)
            return
        end
        tbLocal = forwardTB(hiLocal, hsLocal, logical(isMYI), model);
        if ~isnumeric(tbLocal) || numel(tbLocal) ~= numel(hiLocal) || ...
                any(~isfinite(tbLocal))
            error('invert_pointwise_freeboard_pair:ForwardOutput', ...
                '前向函数必须返回与逐点输入同形的有限 TB。');
        end
        tbLocal = double(tbLocal(:));
        hiMean = sum(weights .* hiLocal);
        hsMean = sum(weights .* hsLocal);
        tbMean = sum(weights .* tbLocal);
        r = tbMean - tbObs;
    end

    function accepted = add_solution(hsValue)
        accepted = false;
        if ~isempty(solutions) && ...
                any(abs([solutions.candidateHs] - hsValue) <= ...
                options.snowToleranceM)
            accepted = true;
            return
        end
        [hiMean, hsMean, tbMean, r, inundatedFraction] = evaluate(hsValue);
        if ~isfinite(r) || abs(r) > options.tbToleranceK, return; end
        solutions(end+1) = struct('hi', hiMean, 'hs', hsMean, ...
            'candidateHs', hsValue, 'tbSim', tbMean, 'tbResidual', r, ...
            'inundatedFraction', inundatedFraction, ...
            'nSamples', numel(fbSamples));
        accepted = true;
    end
end
