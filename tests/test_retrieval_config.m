function tests = test_retrieval_config
%TEST_RETRIEVAL_CONFIG 统一配置入口和派生接口的回归测试。
    tests = functiontests(localfunctions);
end

function testBaselineSchemaAndDerivedNames(testCase)
    config = retrieval_config();
    requiredSections = {'meta','paths','data','collocation','model', ...
        'temperature','domain','inversion','validation','solutionSpace','reference', ...
        'runtime','output'};
    verifyTrue(testCase, all(isfield(config, requiredSections)));
    removedAliases = {'id','years','gridVersion','outputStem', ...
        'brinePolicy','status'};
    verifyFalse(testCase, any(isfield(config, removedAliases)));
    verifyEqual(testCase, config.data.years, [2012 2015]);
    verifyEqual(testCase, config.data.brightnessTemperatureSource, "smos");
    verifyEqual(testCase, config.data.smap.activeProduct, "36km");
    verifyEqual(testCase, config.data.smap.standard.resolutionKm, 36);
    verifyEqual(testCase, config.data.smap.enhanced.resolutionKm, 9);
    verifyEqual(testCase, numel(config.data.smap.oibMatchedDates), 34);
    verifyTrue(testCase, isfolder(config.paths.readSmapRoot));
    verifyEqual(testCase, config.meta.status, "frozen 2026-11-09");
    verifyEqual(testCase, config.meta.frozenAt, "2026-11-09");
    verifyEqual(testCase, config.meta.schemaVersion, 4);
    verifyEqual(testCase, config.data.gridVersion, "EASE2");
    verifyEqual(testCase, config.collocation.windowSize, 3);
    verifyEqual(testCase, config.collocation.centerStride, 3);
    verifyEqual(testCase, config.collocation.centerPhase, [2 2]);
    verifyFalse(testCase, config.model.limitBrinePermittivity);
    verifyEqual(testCase, config.model.radiativeTransferMode, "burke");
    verifyEqual(testCase, config.temperature.source, "fixed");
    verifyEqual(testCase, config.model.surfaceTempC, -30);
    verifyEqual(testCase, config.inversion.maxIceThicknessM, ...
        config.domain.maxIceThicknessM);
    verifyEqual(testCase, config.solutionSpace.maxSnowDepthM, ...
        config.domain.maxSnowDepthM);
    verifyEqual(testCase, config.output.pairsStem, ...
        "pairs_2012_2015_ease2_vb_uncapped_fixed_temp_stride3_same_day");
    verifyEqual(testCase, config.output.inversionStem, ...
        "inversion_2012_2015_ease2_vb_uncapped_fixed_temp_stride3");

    legacy = retrieval_config(struct( ...
        'data', struct('gridVersion', "EASE")));
    verifyTrue(testCase, contains(legacy.output.inversionStem, ...
        "legacy_ease"));
    verifyNotEqual(testCase, legacy.output.inversionStem, ...
        config.output.inversionStem);
end

function testFullConfigCannotBeMixedWithLegacyArguments(testCase)
    config = retrieval_config();
    verifyError(testCase, ...
        @() collocate_day([], [], datetime(2012,3,22), ...
            config, @tb_forward), ...
        'collocate_day:DuplicateConfiguration');
end

function testSmapSourceAndProductOverrides(testCase)
    overrides = struct();
    overrides.data.brightnessTemperatureSource = "smap";
    overrides.data.years = 2016:2019;
    overrides.data.smap.activeProduct = "9km";
    config = retrieval_config(overrides);

    verifyEqual(testCase, config.data.brightnessTemperatureSource, "smap");
    verifyEqual(testCase, config.data.smap.activeProduct, "9km");
    verifyEqual(testCase, config.data.years, 2016:2019);
    verifyEqual(testCase, config.output.inversionStem, ...
        "inversion_2016_2017_2018_2019_" + ...
        "ease2_vb_uncapped_fixed_temp_smap9km_stride3");
end

function testKnownNestedOverridesAndOutputTag(testCase)
    overrides = struct( ...
        'temperature', struct('source', "oib_kt19"), ...
        'domain', struct('maxIceThicknessM', 8));
    config = retrieval_config(overrides);
    verifyFalse(testCase, config.model.limitBrinePermittivity);
    verifyEqual(testCase, config.temperature.source, "oib_kt19");
    verifyEqual(testCase, config.model.brineTemperatureFloorC, -30);
    verifyEqual(testCase, config.inversion.maxIceThicknessM, 8);
    verifyEqual(testCase, config.solutionSpace.maxIceThicknessM, 8);
    verifyEqual(testCase, config.output.inversionStem, ...
        "inversion_2012_2015_ease2_vb_uncapped_oib_kt19_stride3");

    defaults = retrieval_config();
    verifyFalse(testCase, defaults.model.limitBrinePermittivity);
    verifyEqual(testCase, defaults.domain.maxIceThicknessM, 10);

    capConfig = retrieval_config(struct( ...
        'model', struct('limitBrinePermittivity', true)));
    verifyTrue(testCase, capConfig.model.limitBrinePermittivity);
    verifyEqual(testCase, capConfig.output.inversionStem, ...
        "inversion_2012_2015_ease2_vb70_cap_fixed_temp_stride3");
end

function testLinearColdBrineExtrapolationConfig(testCase)
    overrides = struct( ...
        'temperature', struct('source', "oib_kt19", ...
            'oibBrineFloorC', NaN), ...
        'model', struct('brineColdExtrapolationMode', ...
            "linear_below_minus30"));
    config = retrieval_config(overrides);
    verifyTrue(testCase, isnan(config.model.brineTemperatureFloorC));
    verifyEqual(testCase, config.model.brineColdExtrapolationMode, ...
        "linear_below_minus30");
    verifyTrue(testCase, contains(config.output.inversionStem, ...
        "oib_kt19_cold_linear"));
end

function testUnknownOverrideRejected(testCase)
    overrides = struct('model', struct('brineCap', true));
    verifyError(testCase, ...
        @() retrieval_config(overrides), ...
        'retrieval_config:UnknownField');
end

function testInvalidWindowRejected(testCase)
    overrides = struct('collocation', struct('windowSize', 4));
    verifyError(testCase, ...
        @() retrieval_config(overrides), ...
        'retrieval_config:Window');
end

function testSolutionSpaceStructAndLegacyInterfacesAgree(testCase)
    config = retrieval_config();
    coarseOptions = config.solutionSpace;
    coarseOptions.iceStepM = 5;
    coarseOptions.snowStepM = 0.5;

    structSpace = build_solution_space( ...
        false, config.model, @tb_forward, coarseOptions);
    legacySpace = build_solution_space( ...
        false, config.model, @tb_forward, 5, 0.5);

    verifyEqual(testCase, structSpace, legacySpace);
end
