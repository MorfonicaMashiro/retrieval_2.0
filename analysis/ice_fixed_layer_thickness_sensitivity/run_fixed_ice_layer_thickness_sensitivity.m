function result = run_fixed_ice_layer_thickness_sensitivity()
%RUN_FIXED_ICE_LAYER_THICKNESS_SENSITIVITY 测试 1/5/10 cm 固定冰层厚度。
%   本函数只覆盖独立测试分支：不修改 retrieval_config，亦不修改
%   tb_forward 的默认 20 层冰、1 层雪。论文表 1 的冰厚均可被三种目标层厚
%   整除，因此每个场景通过 N=hi/dz 得到严格等厚冰层。

    testRoot = fileparts(mfilename('fullpath'));
    matlabRoot = fileparts(fileparts(testRoot));
    outputRoot = fullfile(testRoot, 'results');
    addpath(matlabRoot);
    if ~isfolder(outputRoot), mkdir(outputRoot); end

    config = retrieval_config();
    reference = config.reference.table1;
    baseModel = config.model;
    baseModel.nSnowLayers = 1;
    baseModel.radiativeTransferMode = "burke";
    baseModel.limitBrinePermittivity = false;

    targetThicknessM = [0.01, 0.05, 0.10];
    targetThicknessCm = 100 .* targetThicknessM;
    scenarioCount = numel(reference.scenario);
    targetCount = numel(targetThicknessM);
    tb20 = nan(scenarioCount, 1);
    for scenarioIndex = 1:scenarioCount
        tb20(scenarioIndex) = tb_forward( ...
            reference.hiM(scenarioIndex), reference.hsM(scenarioIndex), ...
            reference.isMYI(scenarioIndex), baseModel);
    end

    rowCount = scenarioCount * targetCount;
    scenario = strings(rowCount, 1);
    isMYI = false(rowCount, 1);
    hiM = nan(rowCount, 1);
    hsM = nan(rowCount, 1);
    layerThicknessCm = nan(rowCount, 1);
    nIceLayers = nan(rowCount, 1);
    tbPaperK = nan(rowCount, 1);
    tbModelK = nan(rowCount, 1);
    tb20LayerK = nan(rowCount, 1);
    differenceFromPaperK = nan(rowCount, 1);
    differenceFrom20LayerK = nan(rowCount, 1);
    maxBrineVolumePermille = nan(rowCount, 1);
    brineValidityExceeded = false(rowCount, 1);

    rowIndex = 0;
    for targetIndex = 1:targetCount
        dz = targetThicknessM(targetIndex);
        for scenarioIndex = 1:scenarioCount
            rowIndex = rowIndex + 1;
            hi = reference.hiM(scenarioIndex);
            layerCount = round(hi / dz);
            if layerCount < 1 || abs(layerCount * dz - hi) > 1e-10
                error('fixedIceLayerSensitivity:NonDivisibleThickness', ...
                    ['场景 %s 的冰厚 %.6f m 不能被 %.3f m 严格整除；', ...
                    '本测试分支不改变最后一层厚度。'], ...
                    reference.scenario(scenarioIndex), hi, dz);
            end
            model = baseModel;
            model.nIceLayers = layerCount;
            [tb, detail] = tb_forward(hi, reference.hsM(scenarioIndex), ...
                reference.isMYI(scenarioIndex), model);

            scenario(rowIndex) = reference.scenario(scenarioIndex);
            isMYI(rowIndex) = reference.isMYI(scenarioIndex);
            hiM(rowIndex) = hi;
            hsM(rowIndex) = reference.hsM(scenarioIndex);
            layerThicknessCm(rowIndex) = targetThicknessCm(targetIndex);
            nIceLayers(rowIndex) = layerCount;
            tbPaperK(rowIndex) = reference.tbK(scenarioIndex);
            tbModelK(rowIndex) = tb;
            tb20LayerK(rowIndex) = tb20(scenarioIndex);
            differenceFromPaperK(rowIndex) = tb - tbPaperK(rowIndex);
            differenceFrom20LayerK(rowIndex) = tb - tb20(scenarioIndex);
            maxBrineVolumePermille(rowIndex) = ...
                detail.maxBrineVolumePermille;
            brineValidityExceeded(rowIndex) = ...
                detail.brineValidityExceeded;
        end
    end
    details = table(scenario, isMYI, hiM, hsM, layerThicknessCm, ...
        nIceLayers, tbPaperK, tbModelK, tb20LayerK, ...
        differenceFromPaperK, differenceFrom20LayerK, ...
        maxBrineVolumePermille, brineValidityExceeded);

    summaryRows = cell(targetCount, 7);
    for targetIndex = 1:targetCount
        mask = details.layerThicknessCm == targetThicknessCm(targetIndex);
        errorPaper = details.differenceFromPaperK(mask);
        difference20 = details.differenceFrom20LayerK(mask);
        summaryRows(targetIndex, :) = {targetThicknessCm(targetIndex), ...
            min(details.nIceLayers(mask)), max(details.nIceLayers(mask)), ...
            mean(errorPaper), sqrt(mean(errorPaper .^ 2)), ...
            mean(abs(difference20)), max(abs(difference20))};
    end
    summary = cell2table(summaryRows, 'VariableNames', ...
        {'layerThicknessCm','minIceLayers','maxIceLayers', ...
        'biasVsPaperK','rmseVsPaperK','meanAbsDifferenceVs20LayerK', ...
        'maxAbsDifferenceVs20LayerK'});

    wide = table(reference.scenario, reference.isMYI, reference.hiM, ...
        reference.hsM, reference.tbK, tb20, ...
        'VariableNames', {'scenario','isMYI','hiM','hsM','tbPaperK', ...
        'tb20LayerK'});
    for targetIndex = 1:targetCount
        mask = details.layerThicknessCm == targetThicknessCm(targetIndex);
        suffix = sprintf('%02dcm', round(targetThicknessCm(targetIndex)));
        wide.(['nIce_' suffix]) = details.nIceLayers(mask);
        wide.(['tb_' suffix '_K']) = details.tbModelK(mask);
        wide.(['delta_' suffix '_minus20_K']) = ...
            details.differenceFrom20LayerK(mask);
    end

    result = struct('details', details, 'summary', summary, 'wide', wide, ...
        'targetThicknessM', targetThicknessM, 'model', baseModel);
    save(fullfile(outputRoot, ...
        'fixed_ice_layer_thickness_1_5_10cm.mat'), '-struct', 'result');
    writetable(details, fullfile(outputRoot, ...
        'fixed_ice_layer_thickness_1_5_10cm_details.csv'));
    writetable(summary, fullfile(outputRoot, ...
        'fixed_ice_layer_thickness_1_5_10cm_summary.csv'));
    writetable(wide, fullfile(outputRoot, ...
        'fixed_ice_layer_thickness_1_5_10cm_wide.csv'));
    disp(summary);
    fprintf('固定冰层厚度测试结果：%s\n', outputRoot);
end
