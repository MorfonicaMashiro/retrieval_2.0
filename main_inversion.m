% 2017 年论文反演入口。先运行 main 完成同日配准。
% 当前前向辐射模型尚未实现，默认只准备输入；补全后开启 runSolver。

matlabRoot = fileparts(mfilename('fullpath'));
addpath(matlabRoot);
addpath(fullfile(matlabRoot, 'plot'));
cacheFile = fullfile(matlabRoot, 'processed', 'pairs_2017_same_day_v2.mat');
if ~isfile(cacheFile)
    error('main_inversion:MissingPairs', '请先在本目录运行 main 完成配准。');
end
load(cacheFile, 'pairs');

targetYears = 2012; % 可改为 [2012 2015]，仍对每条配对独立反演。
inputs = pairs(ismember(year(pairs.day), targetYears), ...
    {'day','flight','row','col','TB','FBice','FBsnow','isMYI','hi','hs'});
if isempty(inputs), error('main_inversion:NoPairs', '所选年份没有配对。'); end

% 与 2017 论文第 2.1 节默认参数一致；密度为 kg/m^3。
model = struct('surfaceTempC', -30, 'waterTempC', -1.8, ...
    'waterSalinity', 33, 'rhoWater', 1024, ...
    'rhoIce', 915, 'rhoSnow', 320);
fprintf('%s：%d 个有效飞行日、%d 条同日配对。\n', ...
    strjoin(string(targetYears), ','), numel(unique(inputs.day)), height(inputs));
fprintf('每条配对分别使用 TB+FBice 或 TB+FBsnow；OIB hi、hs 留作验证。\n');

runSolver = false; % 完成 tb_forward_2017 的物理亮温计算后改为 true。
if runSolver
    % 先检查前向模型；没有物理亮温函数时明确停止，不输出假解。
    tb_forward_2017(1, 0.1, false, model);
    n = height(inputs);
    hiRadar = nan(n, 1); hsRadar = nan(n, 1);
    hiLaser = nan(n, 1); hsLaser = nan(n, 1);
    hiLaserAlt = nan(n, 1); hsLaserAlt = nan(n, 1);
    nRadarRoots = zeros(n, 1); nLaserRoots = zeros(n, 1);
    tbFromOIB = nan(n, 1);
    for k = 1:n
        [radar, ~] = invert_2017_pair( ...
            inputs.TB(k), inputs.FBice(k), "ice", inputs.isMYI(k), ...
            model, @tb_forward_2017);
        [laser, ~] = invert_2017_pair( ...
            inputs.TB(k), inputs.FBsnow(k), "snow", inputs.isMYI(k), ...
            model, @tb_forward_2017);
        nRadarRoots(k) = numel(radar);
        nLaserRoots(k) = numel(laser);
        if ~isempty(radar)
            hiRadar(k) = radar(1).hi; hsRadar(k) = radar(1).hs;
        end
        if ~isempty(laser)
            % 论文第 5.2 节：双解时默认选雪深较小的一支作验证。
            hiLaser(k) = laser(1).hi; hsLaser(k) = laser(1).hs;
        end
        if numel(laser) > 1
            hiLaserAlt(k) = laser(2).hi; hsLaserAlt(k) = laser(2).hs;
        end
        tbFromOIB(k) = tb_forward_2017( ...
            inputs.hi(k), inputs.hs(k), inputs.isMYI(k), model);
    end
    results = addvars(inputs, hiRadar, hsRadar, hiLaser, hsLaser, ...
        hiLaserAlt, hsLaserAlt, nRadarRoots, nLaserRoots, tbFromOIB);
    fitAll = struct('radarHi', fit_retrieval_2017(results.hi, hiRadar), ...
        'radarHs', fit_retrieval_2017(results.hs, hsRadar), ...
        'laserHi', fit_retrieval_2017(results.hi, hiLaser), ...
        'laserHs', fit_retrieval_2017(results.hs, hsLaser));
    nearTB = abs(tbFromOIB - inputs.TB) <= 1.5;
    fitNear = struct('radarHi', fit_retrieval_2017(results.hi(nearTB), hiRadar(nearTB)), ...
        'radarHs', fit_retrieval_2017(results.hs(nearTB), hsRadar(nearTB)), ...
        'laserHi', fit_retrieval_2017(results.hi(nearTB), hiLaser(nearTB)), ...
        'laserHs', fit_retrieval_2017(results.hs(nearTB), hsLaser(nearTB)));
    spaceFYI = solution_space_2017(false, model, @tb_forward_2017);
    spaceMYI = solution_space_2017(true, model, @tb_forward_2017);
    outputFile = fullfile(matlabRoot, 'processed', 'inversion_2017.mat');
    save(outputFile, 'results', 'fitAll', 'fitNear', 'spaceFYI', 'spaceMYI', 'model', '-v7.3');
    plot_solution_space_2017(spaceFYI, spaceMYI, fullfile(matlabRoot, 'figures'));
    fprintf('已保存反演与解空间：%s\n', outputFile);
else
    fprintf('二分求解器与拟合函数已备好；缺少物理亮温模型，当前停在反演输入准备阶段。\n');
end
