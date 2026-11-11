function tests = test_inversion_closure
%TEST_INVERSION_CLOSURE 八场景端到端回归测试。
    tests = functiontests(localfunctions);
end

function testSelfClosureAndPaperDiagnostics(testCase)
    outputDirectory = fullfile(fileparts(mfilename('fullpath')), 'processed');
    config = retrieval_config();
    result = validate_inversion_closure(outputDirectory);
    self = result.summary(result.summary.mode == "self", :);
    paper = result.summary(result.summary.mode == "paper", :);

    verifyEqual(testCase, height(self), 16);
    verifyEqual(testCase, height(paper), 16);
    % Table 2 的双解数来自论文舍入后的输入，不能强加给本地模型自闭环。
    verifyTrue(testCase, all(self.rootCount >= 1));
    verifyTrue(testCase, all(self.truthRecovered));
    verifyTrue(testCase, all(self.correctBranchOrder));
    verifyLessThanOrEqual(testCase, ...
        max(abs(result.roots.tbResidualK( ...
        result.roots.mode == "self"))), ...
        config.validation.closureTbToleranceK);

    trueBranches = result.roots(result.roots.mode == "self" & ...
        result.roots.rootIndex == 1, :);
    verifyLessThanOrEqual(testCase, ...
        max(abs(trueBranches.hiErrorFromTruthM)), 1e-3);
    verifyLessThanOrEqual(testCase, ...
        max(abs(trueBranches.hsErrorFromTruthM)), 1e-4);
    alternateBranches = result.roots(result.roots.mode == "self" & ...
        result.roots.rootIndex == 2, :);
    verifyTrue(testCase, all(ismember( ...
        alternateBranches.scenario, ...
        config.reference.table2.doubleSolutionScenarios)));
    verifyTrue(testCase, all(alternateBranches.hsErrorFromTruthM > 0));
    verifyTrue(testCase, all(alternateBranches.hiErrorFromTruthM < 0));

    verifyTrue(testCase, all(paper.tbForwardMinusInputK ~= 0));
    noRoot = paper(paper.rootCount == 0, :);
    % 有效解域可能因静力边界分成不连续片段，不能用全局 min/max 判断
    % 必然存在根；无解只要求求解器返回明确状态和零根。
    verifyTrue(testCase, all(noRoot.status == "no_root"));
    verifyTrue(testCase, all(noRoot.rootCount == 0));
end
