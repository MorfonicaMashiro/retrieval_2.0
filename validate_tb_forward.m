function validation = validate_tb_forward(outputDirectory)
%   VALIDATE_TB_FORWARD 对照论文表 1 并检查 1/20/100 层收敛性。
%   VALIDATION = VALIDATE_TB_FORWARD() 返回两个表：
%   validation.scenarios 为 2017 联合反演论文表 1 的 8 个场景；
%   validation.layers 为 1、20、100 层的亮温比较。
%
%   可选 OUTPUTDIRECTORY 非空时，将两张表写为 CSV。这里的误差只用于
%   判断公开公式补全实现与论文数值的接近程度，不进行经验调参。

    if nargin < 1
        outputDirectory = '';
    end
    retrievalConfig = retrieval_config();
    model = retrievalConfig.model;
    % 固定参数来自 Data Synergy 论文 Sect. 2.1（PDF pp. 4--5）：
    % Tsurf=-30 degC, Twater=-1.8 degC, Swater=33 g/kg，密度为
    % rho_water/rho_ice/rho_snow=1024/915/320 kg/m^3。

    table1 = retrievalConfig.reference.table1;
    scenario = table1.scenario;
    isMYI = table1.isMYI;
    hi = table1.hiM;
    hs = table1.hsM;
    tbPaper = table1.tbK;
    % hi、hs、冰型和 TB 均逐项抄录自该论文 Sect. 3、Table 1（PDF p. 10）。
    tbModel = nan(size(hi));
    maxBrinePermille = nan(size(hi));
    brineValidityExceeded = false(size(hi));
    for index = 1:numel(hi)
        [tbModel(index), detail] = tb_forward( ...
            hi(index), hs(index), isMYI(index), model);
        maxBrinePermille(index) = detail.maxBrineVolumePermille;
        brineValidityExceeded(index) = detail.brineValidityExceeded;
    end
    differenceK = tbModel - tbPaper;
    scenarios = table(scenario, isMYI, hi, hs, tbPaper, tbModel, ...
        differenceK, maxBrinePermille, brineValidityExceeded);

    layerSensitivity = retrievalConfig.validation.layerSensitivity;
    testHi = layerSensitivity.hiM;
    testHs = layerSensitivity.hsM;
    testMYI = layerSensitivity.isMYI;
    layerCounts = layerSensitivity.iceLayerCounts;
    % 层数敏感性设置来自用户提供的 Improving L-band radiation model 论文
    % Sect. 3（PDF p. 5）：N=1 与 20 差异大，N=20 与 100 的 TB 差异<0.1 K。
    layerResults = cell(size(layerCounts));
    for layerIndex = 1:numel(layerCounts)
        layerModel = model;
        layerModel.nIceLayers = layerCounts(layerIndex);
        current = nan(size(testHi));
        for index = 1:numel(testHi)
            current(index) = tb_forward( ...
                testHi(index), testHs(index), testMYI(index), layerModel);
        end
        layerResults{layerIndex} = current;
    end
    tb1 = layerResults{1};
    tb20 = layerResults{2};
    tb100 = layerResults{3};
    layers = table(testMYI, testHi, testHs, tb1, tb20, tb100, ...
        tb20 - tb1, tb100 - tb20, ...
        'VariableNames', {'isMYI','hi','hs','TB_1layer','TB_20layers', ...
        'TB_100layers','delta20minus1','delta100minus20'});

    validation = struct('scenarios', scenarios, 'layers', layers, ...
        'scenarioRMSEK', sqrt(mean(differenceK.^2)), ...
        'scenarioBiasK', mean(differenceK), ...
        'maxAbs20to100K', max(abs(layers.delta100minus20)), ...
        'retrievalConfig', retrievalConfig);

    if strlength(string(outputDirectory)) > 0
        if ~isfolder(outputDirectory), mkdir(outputDirectory); end
        outputStem = retrievalConfig.output.forwardValidationStem;
        writetable(scenarios, fullfile(outputDirectory, ...
            outputStem + "_table1.csv"));
        writetable(layers, fullfile(outputDirectory, ...
            outputStem + "_layers.csv"));
    end
end
