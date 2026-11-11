function space = build_solution_space( ...
        isMYI, configOrModel, forwardTB, optionsOrStepIceM, legacyStepSnowM)
%BUILD_SOLUTION_SPACE 扫描冰厚、雪深、干舷和亮温解空间。
%   SPACE = BUILD_SOLUTION_SPACE(ISMYI,CONFIG) 是推荐接口。
%   旧式 MODEL、FORWARDTB、OPTIONS/步长输入仍兼容。
%   第四个输入推荐传入 config.solutionSpace 结构。为兼容旧脚本，第四、
%   第五个输入仍可分别传数值 iceStepM、snowStepM。

    config = retrieval_config();
    isConfigInterface = nargin >= 2 && isstruct(configOrModel) && ...
        all(isfield(configOrModel, {'model','solutionSpace'}));
    if isConfigInterface
        config = configOrModel;
        model = config.model;
        forwardTB = @tb_forward;
        options = config.solutionSpace;
    else
        model = configOrModel;
        if nargin < 3 || isempty(forwardTB), forwardTB = @tb_forward; end
        options = config.solutionSpace;
    end
    if ~isConfigInterface && nargin >= 4 && ~isempty(optionsOrStepIceM)
        if isstruct(optionsOrStepIceM)
            options = apply_known_options(options, optionsOrStepIceM);
        elseif isnumeric(optionsOrStepIceM) && isscalar(optionsOrStepIceM)
            options.iceStepM = optionsOrStepIceM;
            if nargin >= 5 && ~isempty(legacyStepSnowM)
                options.snowStepM = legacyStepSnowM;
            end
        else
            error('build_solution_space:Options', ...
                '第四个输入必须是配置结构或冰厚扫描步长。');
        end
    end
    validate_inputs(isMYI, model, forwardTB, options);

    iceAxis = options.minIceThicknessM:options.iceStepM: ...
        options.maxIceThicknessM;
    snowAxis = options.minSnowDepthM:options.snowStepM: ...
        options.maxSnowDepthM;
    [hiGrid, hsGrid] = ndgrid(iceAxis, snowAxis);
    rhoW = model.rhoWater;
    rhoI = model.rhoIce;
    rhoS = model.rhoSnow;
    fbIce = ((rhoW - rhoI) .* hiGrid - rhoS .* hsGrid) ./ rhoW;
    fbSnow = ((rhoW - rhoI) .* hiGrid + (rhoW - rhoS) .* hsGrid) ./ rhoW;

    valid = true(size(fbIce));
    if options.excludeFloodedIce
        % 论文第 2.2 节排除 FBice<0 的淹没冰情形。
        valid = valid & fbIce >= 0;
    end
    hi = hiGrid(valid);
    hs = hsGrid(valid);
    fbIce = fbIce(valid);
    fbSnow = fbSnow(valid);
    tb = nan(size(hi));
    for index = 1:numel(hi)
        tb(index) = forwardTB(hi(index), hs(index), isMYI, model);
    end
    if any(~isfinite(tb))
        error('build_solution_space:ForwardOutput', ...
            '前向模型返回无效亮温。');
    end
    space = table(hi, hs, fbIce, fbSnow, tb, ...
        repmat(logical(isMYI), numel(hi), 1), ...
        'VariableNames', {'hi','hs','FBice','FBsnow','TB','isMYI'});
end

function options = apply_known_options(options, overrides)
    names = fieldnames(overrides);
    unknown = setdiff(names, fieldnames(options));
    if ~isempty(unknown)
        error('build_solution_space:Options', ...
            '未知解空间选项：%s。', strjoin(unknown, ', '));
    end
    for index = 1:numel(names)
        options.(names{index}) = overrides.(names{index});
    end
end

function validate_inputs(isMYI, model, forwardTB, options)
    requiredModel = {'rhoWater','rhoIce','rhoSnow'};
    numericOptions = {'minIceThicknessM','maxIceThicknessM','iceStepM', ...
        'minSnowDepthM','maxSnowDepthM','snowStepM'};
    for index = 1:numel(numericOptions)
        value = options.(numericOptions{index});
        if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value)
            error('build_solution_space:Inputs', ...
                '解空间数值参数必须是有限标量。');
        end
    end
    if ~isscalar(isMYI) || ~ismember(double(isMYI), [0 1]) || ...
            ~isstruct(model) || ~all(isfield(model, requiredModel)) || ...
            ~isa(forwardTB, 'function_handle') || ...
            options.minIceThicknessM < 0 || ...
            options.maxIceThicknessM <= options.minIceThicknessM || ...
            options.iceStepM <= 0 || options.minSnowDepthM < 0 || ...
            options.maxSnowDepthM <= options.minSnowDepthM || ...
            options.snowStepM <= 0 || ...
            ~(isscalar(options.excludeFloodedIce) && ...
            (islogical(options.excludeFloodedIce) || ...
            ismember(options.excludeFloodedIce, [0 1])))
        error('build_solution_space:Inputs', ...
            '冰型、前向模型或解空间配置无效。');
    end
end
