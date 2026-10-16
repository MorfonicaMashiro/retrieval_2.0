function [solutions, info] = invert_2017_pair(tbObs, fbObs, branch, isMYI, model, forwardTB, options)
%INVERT_2017_PAIR 按 2017 论文图 5/6 求一条配对的所有物理解。
%   branch="ice"：TB + 冰面干舷，先以 4 cm 步长寻找变号区间，再二分。
%   branch="snow"：TB + 雪面干舷，先在 0~1 m 以 1 cm 扫描，再逐区间二分。
%   forwardTB(hi,hs,isMYI,model) 必须返回模拟亮温，单位 K。
%   solutions 为结构数组；雪面干舷路线可能返回两个解。
%   未找到解时 solutions 为空，info.status="no_root"。

    if nargin < 7 || isempty(options), options = struct(); end
    defaults = struct('maxSnowDepthM', 1, 'maxIceThicknessM', 10, ...
        'initialStepM', 0.04, 'scanStepM', 0.01, ...
        'tbToleranceK', 0.01, 'snowToleranceM', 1e-5, ...
        'maxIterations', 60);
    names = fieldnames(defaults);
    for j = 1:numel(names)
        if ~isfield(options, names{j})
            options.(names{j}) = defaults.(names{j});
        end
    end

    branch = string(branch);
    required = {'rhoWater', 'rhoIce', 'rhoSnow'};
    if ~isscalar(tbObs) || ~isfinite(tbObs) || ...
            ~isscalar(fbObs) || ~isfinite(fbObs) || ...
            ~ismember(branch, ["ice", "snow"]) || ...
            ~isscalar(isMYI) || ~ismember(double(isMYI), [0 1]) || ...
            ~isstruct(model) || ~all(isfield(model, required)) || ...
            ~isa(forwardTB, 'function_handle')
        error('invert_2017_pair:Inputs', '观测、冰型、密度或前向函数无效。');
    end
    rhoW = model.rhoWater;
    rhoI = model.rhoIce;
    rhoS = model.rhoSnow;
    if ~(rhoW > rhoI && rhoI > rhoS && rhoS > 0) || ...
            options.maxSnowDepthM <= 0 || options.maxIceThicknessM <= 0 || ...
            options.initialStepM <= 0 || options.scanStepM <= 0 || ...
            options.tbToleranceK <= 0 || options.snowToleranceM <= 0 || ...
            options.maxIterations < 1
        error('invert_2017_pair:Parameters', '密度或求解配置不合理。');
    end

    solutions = struct('hi', {}, 'hs', {}, 'tbSim', {}, 'tbResidual', {});
    info = struct('status', "no_root", 'scanCount', 0, ...
        'rootCount', 0, 'unconvergedBrackets', 0);
    if fbObs < 0
        info.status = "negative_freeboard";
        return
    end

    hsLimit = options.maxSnowDepthM;
    if branch == "snow"
        % 论文第 2.2 节排除淹没海冰：FBice = FBsnow - hs >= 0。
        hsLimit = min(hsLimit, fbObs);
        step = options.scanStepM;
    else
        step = options.initialStepM;
    end
    hsGrid = 0:step:hsLimit;
    if isempty(hsGrid) || hsGrid(end) < hsLimit
        hsGrid(end+1) = hsLimit;
    end
    residual = nan(size(hsGrid));
    for j = 1:numel(hsGrid)
        [~, ~, residual(j)] = evaluate(hsGrid(j));
    end
    info.scanCount = numel(hsGrid);

    % 处理扫描点恰好命中阈值的情形；连续命中视为一个平坦区间。
    hit = isfinite(residual) & abs(residual) <= options.tbToleranceK;
    edges = diff([false, hit, false]);
    first = find(edges == 1);
    last = find(edges == -1) - 1;
    for j = 1:numel(first)
        add_solution(hsGrid(round((first(j) + last(j)) / 2)));
    end

    % 相邻模拟 TB 分居观测 TB 两侧时，按图 5/6 缩小区间。
    for j = 1:numel(hsGrid)-1
        if hit(j) || hit(j+1) || ...
                ~isfinite(residual(j)) || ~isfinite(residual(j+1)) || ...
                residual(j) * residual(j+1) >= 0
            continue
        end
        lo = hsGrid(j); hiBracket = hsGrid(j+1);
        rLo = residual(j);
        converged = false;
        for iter = 1:options.maxIterations
            middle = (lo + hiBracket) / 2;
            [~, ~, rMiddle] = evaluate(middle);
            if abs(rMiddle) <= options.tbToleranceK
                converged = add_solution(middle);
                break
            end
            if rLo * rMiddle <= 0
                hiBracket = middle;
            else
                lo = middle;
                rLo = rMiddle;
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
        [~, order] = sort([solutions.hs]);
        solutions = solutions(order);
        info.rootCount = numel(solutions);
        info.status = "solved";
        if branch == "ice" && info.rootCount > 1
            info.status = "multiple_radar_roots";
        end
    end

    function [hiValue, tbValue, r] = evaluate(hsValue)
        if branch == "ice"
            % 论文式 (2)。
            hiValue = (rhoW * fbObs + rhoS * hsValue) / (rhoW - rhoI);
            fbIce = fbObs;
        else
            % 论文式 (3)。
            hiValue = (rhoW * fbObs - (rhoW - rhoS) * hsValue) ...
                / (rhoW - rhoI);
            fbIce = fbObs - hsValue;
        end
        tbValue = NaN; r = NaN;
        if hiValue < 0 || hiValue > options.maxIceThicknessM || fbIce < 0
            return
        end
        tbValue = forwardTB(hiValue, hsValue, isMYI, model);
        if ~isnumeric(tbValue) || ~isscalar(tbValue) || ~isfinite(tbValue)
            error('invert_2017_pair:ForwardOutput', ...
                '前向函数必须返回有限的标量亮温。');
        end
        r = tbValue - tbObs;
    end

    function accepted = add_solution(hsValue)
        accepted = false;
        if ~isempty(solutions) && ...
                any(abs([solutions.hs] - hsValue) <= options.snowToleranceM)
            accepted = true;
            return
        end
        [hiValue, tbValue, r] = evaluate(hsValue);
        if ~isfinite(r) || abs(r) > options.tbToleranceK
            return
        end
        solutions(end+1) = struct('hi', hiValue, 'hs', hsValue, ...
            'tbSim', tbValue, 'tbResidual', r); %#ok<AGROW>
        accepted = true;
    end
end
