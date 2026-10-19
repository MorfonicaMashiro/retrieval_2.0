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
    model = struct('surfaceTempC', -30, 'waterTempC', -1.8, ...
        'waterSalinity', 33, 'rhoWater', 1024, ...
        'rhoIce', 915, 'rhoSnow', 320);
    % 固定参数来自 Data Synergy 论文 Sect. 2.1（PDF pp. 4--5）：
    % Tsurf=-30 degC, Twater=-1.8 degC, Swater=33 g/kg，密度为
    % rho_water/rho_ice/rho_snow=1024/915/320 kg/m^3。

    scenario = ["I"; "II"; "III"; "IV"; "V"; "VI"; "VII"; "VIII"];
    isMYI = [false; false; false; false; true; true; true; true];
    hi = [0.5; 1.0; 1.5; 2.5; 1.5; 2.5; 3.0; 5.0];
    hs = [0.05; 0.03; 0.10; 0.25; 0.15; 0.15; 0.35; 0.40];
    tbPaper = [231.5; 238.2; 243.7; 246.5; 243.8; 248.9; 252.0; 250.4];
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

    testHi = [0.5; 1.0; 2.5; 5.0; 0.5; 1.5; 3.0; 5.0];
    testHs = [0.05; 0.10; 0.25; 0.40; 0.05; 0.15; 0.35; 0.40];
    testMYI = [false; false; false; false; true; true; true; true];
    layerCounts = [1, 20, 100];
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
        'maxAbs20to100K', max(abs(layers.delta100minus20)));

    if strlength(string(outputDirectory)) > 0
        if ~isfolder(outputDirectory), mkdir(outputDirectory); end
        writetable(scenarios, fullfile(outputDirectory, ...
            'tb_forward_table1_validation.csv'));
        writetable(layers, fullfile(outputDirectory, ...
            'tb_forward_layer_sensitivity.csv'));
    end
end
