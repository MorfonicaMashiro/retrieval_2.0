function [operatorMetrics, retrievalResults] = ...
        summarize_retrieval_with_sample_tb_mean()
%SUMMARIZE_RETRIEVAL_WITH_SAMPLE_TB_MEAN 用两种 TB 算子重算验证子集。
%   反演根本身由 SMOS TB 与平均干舷决定，不随 OIB 验证 TB 的计算次序
%   改变。本函数复用历史求解根，只重算 |TBmodel-TBSMOS|<=1.5 K 掩膜。

    operatorRoot = fileparts(mfilename('fullpath'));
    matlabRoot = fileparts(fileparts(fileparts(operatorRoot)));
    resultRoot = fullfile(operatorRoot, 'results');
    addpath(matlabRoot, operatorRoot);

    inversionFile = fullfile(matlabRoot, 'processed', ...
        'inversion_2012_2015.mat');
    averagingFile = fullfile(resultRoot, ...
        'tb_sample_mean_2012_2015.mat');
    if ~isfile(inversionFile) || ~isfile(averagingFile)
        error('summarize_retrieval_with_sample_tb_mean:MissingInput', ...
            '缺少历史反演结果或本轮 TB 平均结果。');
    end
    inversion = load(inversionFile, 'results');
    averaging = load(averagingFile, 'pairs');
    retrievalResults = inversion.results( ...
        ismember(year(inversion.results.day), [2012 2015]), :);
    pairs = averaging.pairs;
    retrievalResults = sortrows(retrievalResults, {'day','row','col'});
    pairs = sortrows(pairs, {'day','row','col'});
    if height(retrievalResults) ~= height(pairs) || ...
            ~isequal(retrievalResults.day, pairs.day) || ...
            ~isequal(retrievalResults.row, pairs.row) || ...
            ~isequal(retrievalResults.col, pairs.col)
        error('summarize_retrieval_with_sample_tb_mean:PairKeys', ...
            '历史反演与本轮配对的 day/row/col 不一致。');
    end
    inputFields = {'TB','FBice','FBsnow','isMYI','hi','hs','nOIB'};
    for fieldIndex = 1:numel(inputFields)
        field = inputFields{fieldIndex};
        if ~isequaln(retrievalResults.(field), pairs.(field))
            error('summarize_retrieval_with_sample_tb_mean:PairInputs', ...
                '历史反演与本轮配对的 %s 不一致。', field);
        end
    end

    retrievalResults.TBOIBMeanParameters = pairs.TBOIBMeanParameters;
    retrievalResults.TBOIBMeanSamples = pairs.TBOIBMeanSamples;
    retrievalResults.TBSampleMinusParameter = ...
        pairs.TBSampleMinusParameter;
    operatorNames = ["mean_parameters", "mean_samples"];
    operatorFields = {'TBOIBMeanParameters', 'TBOIBMeanSamples'};
    operatorMetrics = table();
    allRows = true(height(retrievalResults), 1);
    for operatorIndex = 1:numel(operatorNames)
        current = retrievalResults;
        current.tbFromOIB = current.(operatorFields{operatorIndex});
        nearTB = abs(current.tbFromOIB - current.TB) <= 1.5;
        metrics = summarize_retrieval_subset(current, allRows, nearTB);
        metrics = addvars(metrics, ...
            repmat(operatorNames(operatorIndex), height(metrics), 1), ...
            nnz(nearTB) .* ones(height(metrics), 1), ...
            'Before', 1, 'NewVariableNames', {'operator','nearTBRows'});
        operatorMetrics = [operatorMetrics; metrics]; %#ok<AGROW>
    end

    outputStem = fullfile(resultRoot, ...
        'retrieval_2012_2015_tb_operator_validation');
    save(outputStem + ".mat", 'operatorMetrics', 'retrievalResults', '-v7.3');
    writetable(operatorMetrics, outputStem + "_metrics.csv");
    fprintf('已用两种 TB 算子重算 2012/2015 反演验证子集。\n');
end
