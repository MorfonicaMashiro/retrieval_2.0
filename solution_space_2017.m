function space = solution_space_2017(isMYI, model, forwardTB, stepIceM, stepSnowM)
%   SOLUTION_SPACE_2017 按论文第 2.2 节扫描海冰厚度、雪深和两种干舷。
%   亮温必须由 forwardTB(hi,hs,isMYI,model) 提供；不使用经验替代值。
%   原论文只给出扫描范围，未规定扫描间隔；默认冰厚 5 cm、雪深 1 cm。

    if nargin < 4 || isempty(stepIceM), stepIceM = 0.05; end
    if nargin < 5 || isempty(stepSnowM), stepSnowM = 0.01; end
    if ~isscalar(isMYI) || ~ismember(double(isMYI), [0 1]) || ...
            ~isstruct(model) || ...
            ~all(isfield(model, {'rhoWater','rhoIce','rhoSnow'})) || ...
            ~isa(forwardTB, 'function_handle') || ...
            ~isscalar(stepIceM) || ~isscalar(stepSnowM) || ...
            ~isfinite(stepIceM) || ~isfinite(stepSnowM) || ...
            stepIceM <= 0 || stepSnowM <= 0
        error('solution_space_2017:Inputs', '冰型、前向模型或扫描间隔无效。');
    end

    [hiGrid, hsGrid] = ndgrid(0:stepIceM:10, 0:stepSnowM:1);
    rhoW = model.rhoWater;
    rhoI = model.rhoIce;
    rhoS = model.rhoSnow;
    fbIce = ((rhoW - rhoI) .* hiGrid - rhoS .* hsGrid) ./ rhoW;
    fbSnow = ((rhoW - rhoI) .* hiGrid + (rhoW - rhoS) .* hsGrid) ./ rhoW;

    % 论文第 2.2 节排除 FBice<0 的淹没冰情形。
    valid = fbIce >= 0;
    hi = hiGrid(valid);
    hs = hsGrid(valid);
    fbIce = fbIce(valid);
    fbSnow = fbSnow(valid);
    tb = nan(size(hi));
    for k = 1:numel(hi)
        tb(k) = forwardTB(hi(k), hs(k), isMYI, model);
    end
    if any(~isfinite(tb))
        error('solution_space_2017:ForwardOutput', '前向模型返回无效亮温。');
    end
    space = table(hi, hs, fbIce, fbSnow, tb, ...
        repmat(logical(isMYI), numel(hi), 1), ...
        'VariableNames', {'hi','hs','FBice','FBsnow','TB','isMYI'});
end
