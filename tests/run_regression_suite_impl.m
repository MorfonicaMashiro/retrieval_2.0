function report = run_regression_suite_impl()
%RUN_REGRESSION_SUITE_IMPL 执行当前冻结版本的完整 MATLAB 回归测试。
%   REPORT = RUN_REGRESSION_SUITE_IMPL() 检查冻结配置，运行单元测试、
%   四条空间平均/反演算子验证，并要求全部 MATLAB 文件通过静态检查。

    testRoot = fileparts(mfilename('fullpath'));
    matlabRoot = fileparts(testRoot);
    addpath(genpath(matlabRoot));

    config = retrieval_config();
    verify_frozen_defaults(config);

    testListing = dir(fullfile(testRoot, 'test_*.m'));
    [~, order] = sort({testListing.name});
    testListing = testListing(order);
    testFiles = fullfile({testListing.folder}, {testListing.name});
    unitResults = runtests(testFiles);
    if ~all([unitResults.Passed])
        disp(table(unitResults));
        error('run_regression_suite_impl:UnitTests', ...
            '单元测试未全部通过，冻结验收中止。');
    end

    verifierNames = [
        "verify_stride3_parameter_mean"
        "verify_sample_tb_mean_operator"
        "verify_weighted_parameter_mean_operator"
        "verify_pointwise_freeboard_inversion"
    ];
    for index = 1:numel(verifierNames)
        feval(char(verifierNames(index)));
    end

    matlabFiles = dir(fullfile(matlabRoot, '**', '*.m'));
    staticMessageCount = 0;
    for index = 1:numel(matlabFiles)
        filePath = fullfile(matlabFiles(index).folder, matlabFiles(index).name);
        messages = checkcode(filePath, '-id');
        if ~isempty(messages)
            staticMessageCount = staticMessageCount + numel(messages);
            fprintf(2, '静态检查消息：%s\n', filePath);
            for messageIndex = 1:numel(messages)
                fprintf(2, '  L%d [%s] %s\n', ...
                    messages(messageIndex).line, messages(messageIndex).id, ...
                    messages(messageIndex).message);
            end
        end
    end
    if staticMessageCount > 0
        error('run_regression_suite_impl:StaticCheck', ...
            'MATLAB 静态检查发现 %d 条消息，冻结验收中止。', ...
            staticMessageCount);
    end

    report = struct( ...
        'passed', true, ...
        'frozenAt', config.meta.frozenAt, ...
        'configurationId', config.meta.id, ...
        'unitTestCount', numel(unitResults), ...
        'verifierCount', numel(verifierNames), ...
        'staticFileCount', numel(matlabFiles), ...
        'staticMessageCount', staticMessageCount, ...
        'pairsStem', config.output.pairsStem, ...
        'inversionStem', config.output.inversionStem);

    fprintf(['\n冻结版本回归测试通过：%d 个单元测试、%d 个算子验证、' ...
        '%d 个 MATLAB 文件静态检查（0 条消息）。\n'], ...
        report.unitTestCount, report.verifierCount, report.staticFileCount);
    fprintf('默认配置：%s\n', report.configurationId);
    fprintf('默认反演输出：%s\n', report.inversionStem);
end

function verify_frozen_defaults(config)
    assert(config.collocation.windowSize == 3);
    assert(config.collocation.centerStride == 3);
    assert(isequal(config.collocation.centerPhase, [2 2]));
    assert(config.model.radiativeTransferMode == "burke");
    assert(config.model.nIceLayers == 20);
    assert(config.model.nSnowLayers == 1);
    assert(~config.model.limitBrinePermittivity);
    assert(config.temperature.source == "fixed");
    assert(config.temperature.fixedC == -30);
    assert(config.temperature.oibBrineFloorC == -30);
    assert(config.data.gridVersion == "EASE2");
end
