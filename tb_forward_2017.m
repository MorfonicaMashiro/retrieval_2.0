function tb = tb_forward_2017(hi, hs, isMYI, model)
%TB_FORWARD_2017 2017 论文式 (1) 的 L 波段亮温模型接口。
%   tb = tb_forward_2017(hi, hs, isMYI, model)
%   hi、hs 的单位为 m；isMYI 为 0(FYI) 或 1(MYI)；输出 tb 的单位为 K。
%   输出应与 SMOS L3B 观测亮温的层级及极化/角度平均约定一致。
%   model.surfaceTempC、waterTempC 为摄氏度，waterSalinity 为 g/kg。
%   论文 [32] 的多层辐射传输、温盐剖面等尚未实现，故不会返回假亮温。

    required = {'surfaceTempC', 'waterTempC', 'waterSalinity'};
    if ~isstruct(model) || ~all(isfield(model, required))
        error('tb_forward_2017:Parameters', ...
            'model 需包含 surfaceTempC、waterTempC、waterSalinity。');
    end
    if ~isequal(size(hi), size(hs)) || ...
            ~isnumeric(hi) || ~isnumeric(hs) || ...
            any(~isfinite(hi(:)) | hi(:) < 0) || ...
            any(~isfinite(hs(:)) | hs(:) < 0) || ...
            ~isscalar(isMYI) || ~ismember(double(isMYI), [0 1])
        error('tb_forward_2017:Inputs', ...
            'hi、hs 需为同形状非负有限数值；isMYI 需为 0 或 1。');
    end

    % Zhou 等 (2017) 文献 [32] 已给出分层、盐度和混合面积比例等约束。
    % 当前还缺少完整的冰面辐射亮温数值实现，不可用混合公式替代。
    % 2017 联合反演论文只给出 TB(hi,hs,Tsfc,ice_type,Twater,Swater)，
    % 不能从这一符号式直接计算数值亮温。
    error('tb_forward_2017:NotImplemented', ...
        'L 波段辐射模型尚未实现；请先补全论文 [32] 的公式与参数。');
end
