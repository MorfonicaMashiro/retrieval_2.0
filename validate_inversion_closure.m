function validation = validate_inversion_closure(outputDirectory, modelOverrides)
%VALIDATE_INVERSION_CLOSURE 八个理想场景的双路线端到端闭环。
%   自闭环使用精确静力平衡干舷和本地 tb_forward 亮温；论文闭环逐字使用
%   Xu et al. (2017) Table 1 的四位小数干舷与一位小数亮温。
%   两组结果分别保存，绝不把论文亮温与本地前向模型的差异算作求解误差。
%
%   VALIDATION = VALIDATE_INVERSION_CLOSURE() 返回 summary、roots
%   两张表。指定 OUTPUTDIRECTORY 时写入同名 CSV。此函数只求解 8 个
%   理想场景，不读取 OIB/SMOS 数据，也不运行参数或层数敏感性分析。

    if nargin < 1, outputDirectory = ''; end
    if nargin < 2, modelOverrides = struct(); end
    retrievalConfig = retrieval_config( ...
        struct('model', modelOverrides));
    model = retrievalConfig.model;
    tolerance = retrievalConfig.validation;
    table1 = retrievalConfig.reference.table1;
    table2 = retrievalConfig.reference.table2;
    scenario = table1.scenario;
    isMYI = table1.isMYI;
    hiTrue = table1.hiM;
    hsTrue = table1.hsM;
    fbIcePaper = table1.fbIceM;
    fbSnowPaper = table1.fbSnowM;
    tbPaper = table1.tbK;
    alternateHi = table2.alternateHiM;
    alternateHs = table2.alternateHsM;

    tbSelf = nan(8, 1);
    for k = 1:8
        tbSelf(k) = tb_forward(hiTrue(k), hsTrue(k), isMYI(k), model);
    end
    fbIceSelf = ((model.rhoWater - model.rhoIce) .* hiTrue ...
        - model.rhoSnow .* hsTrue) ./ model.rhoWater;
    fbSnowSelf = fbIceSelf + hsTrue;

    summaryRows = cell(0, 21);
    rootRows = cell(0, 18);
    routes = ["ice", "snow"];
    modes = ["self", "paper"];
    for mode = modes
        for k = 1:8
            for route = routes
                if mode == "self"
                    inputTB = tbSelf(k);
                    if route == "ice", inputFB = fbIceSelf(k);
                    else, inputFB = fbSnowSelf(k); end
                else
                    inputTB = tbPaper(k);
                    if route == "ice", inputFB = fbIcePaper(k);
                    else, inputFB = fbSnowPaper(k); end
                end
                [solutions, info] = invert_pair(inputTB, inputFB, ...
                    route, isMYI(k), retrievalConfig);
                nearestHiError = NaN;
                nearestHsError = NaN; nearestTBResidual = NaN;
                if ~isempty(solutions)
                    distances = hypot([solutions.hi] - hiTrue(k), ...
                        [solutions.hs] - hsTrue(k));
                    [~, nearest] = min(distances);
                    nearestHiError = solutions(nearest).hi - hiTrue(k);
                    nearestHsError = solutions(nearest).hs - hsTrue(k);
                    nearestTBResidual = solutions(nearest).tbResidual;
                end
                expectedRoots = 1;
                if route == "snow" && ismember( ...
                        scenario(k), table2.doubleSolutionScenarios)
                    expectedRoots = 2;
                end
                hasCorrectOrder = true;
                if numel(solutions) > 1
                    hasCorrectOrder = all(diff([solutions.hs]) > 0);
                    if route == "snow"
                        hasCorrectOrder = hasCorrectOrder && ...
                            all(diff([solutions.hi]) < 0);
                    end
                end
                truthRecovered = ~isempty(solutions) && ...
                    abs(nearestHiError) <= tolerance.closureHiToleranceM && ...
                    abs(nearestHsError) <= tolerance.closureHsToleranceM && ...
                    abs(nearestTBResidual) <= tolerance.closureTbToleranceK;
                summaryRows(end+1, :) = {mode, scenario(k), route, ...
                    isMYI(k), hiTrue(k), hsTrue(k), inputTB, inputFB, ...
                    tbSelf(k), tbSelf(k) - inputTB, expectedRoots, ...
                    info.rootCount, info.status, info.scanCount, ...
                    info.scanTBMinK, info.scanTBMaxK, ...
                    nearestHiError, nearestHsError, nearestTBResidual, ...
                    truthRecovered, hasCorrectOrder}; %#ok<AGROW>

                for j = 1:numel(solutions)
                    expectedHi = NaN; expectedHs = NaN;
                    if mode == "paper"
                        if j == 1
                            expectedHi = hiTrue(k);
                            expectedHs = hsTrue(k);
                        elseif route == "snow" && j == 2
                            expectedHi = alternateHi(k);
                            expectedHs = alternateHs(k);
                        end
                    end
                    rootRows(end+1, :) = {mode, scenario(k), route, ...
                        j, info.rootCount, info.status, inputTB, inputFB, ...
                        solutions(j).hi, solutions(j).hs, ...
                        solutions(j).tbSim, solutions(j).tbResidual, ...
                        solutions(j).hi - hiTrue(k), ...
                        solutions(j).hs - hsTrue(k), expectedHi, ...
                        expectedHs, solutions(j).hi - expectedHi, ...
                        solutions(j).hs - expectedHs}; %#ok<AGROW>
                end
            end
        end
    end
    summary = cell2table(summaryRows, 'VariableNames', ...
        {'mode','scenario','route','isMYI','hiTrueM','hsTrueM', ...
        'tbInputK','freeboardInputM','tbForwardAtTruthK', ...
        'tbForwardMinusInputK','expectedRootCount','rootCount', ...
        'status','scanCount','scanTBMinK','scanTBMaxK', ...
        'nearestHiErrorM','nearestHsErrorM', ...
        'nearestTBResidualK','truthRecovered','correctBranchOrder'});
    roots = cell2table(rootRows, 'VariableNames', ...
        {'mode','scenario','route','rootIndex','rootCount','status', ...
        'tbInputK','freeboardInputM','hiM','hsM','tbSimK', ...
        'tbResidualK','hiErrorFromTruthM','hsErrorFromTruthM', ...
        'paperTable2HiM','paperTable2HsM','deltaFromTable2HiM', ...
        'deltaFromTable2HsM'});
    validation = struct('summary', summary, 'roots', roots, ...
        'model', model, 'retrievalConfig', retrievalConfig);

    if strlength(string(outputDirectory)) > 0
        if ~isfolder(outputDirectory), mkdir(outputDirectory); end
        outputStem = retrievalConfig.output.closureStem;
        writetable(summary, fullfile(outputDirectory, ...
            outputStem + "_summary.csv"));
        writetable(roots, fullfile(outputDirectory, ...
            outputStem + "_roots.csv"));
    end
end
