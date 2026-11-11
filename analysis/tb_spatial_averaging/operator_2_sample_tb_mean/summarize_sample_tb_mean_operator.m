function [metrics, differences] = summarize_sample_tb_mean_operator(pairs)
%SUMMARIZE_SAMPLE_TB_MEAN_OPERATOR 汇总逐样本 TB 平均算子的指标。

    required = {'day','TB','nOIB','TBOIBMeanParameters', ...
        'TBOIBMeanSamples','TBSampleMinusParameter'};
    if ~all(ismember(required, pairs.Properties.VariableNames))
        error('summarize_sample_tb_mean_operator:Fields', ...
            '输入配对表缺少 TB 平均试验字段。');
    end

    groupNames = "all";
    groupMasks = {true(height(pairs), 1)};
    years = unique(year(pairs.day)).';
    for currentYear = years
        groupNames(end + 1, 1) = "year_" + string(currentYear); %#ok<AGROW>
        groupMasks{end + 1, 1} = year(pairs.day) == currentYear; %#ok<AGROW>
    end
    days = unique(pairs.day).';
    for currentDay = days
        groupNames(end + 1, 1) = "day_" + ...
            string(currentDay, 'yyyyMMdd'); %#ok<AGROW>
        groupMasks{end + 1, 1} = pairs.day == currentDay; %#ok<AGROW>
    end

    subsetNames = ["all", "nOIB_ge_50"];
    operatorNames = ["mean_parameters", "mean_samples"];
    operatorFields = {'TBOIBMeanParameters', 'TBOIBMeanSamples'};
    metricRows = cell(0, 11);
    differenceRows = cell(0, 12);
    for groupIndex = 1:numel(groupNames)
        for subsetIndex = 1:numel(subsetNames)
            mask = groupMasks{groupIndex};
            if subsetNames(subsetIndex) == "nOIB_ge_50"
                mask = mask & pairs.nOIB >= 50;
            end
            for operatorIndex = 1:numel(operatorNames)
                fit = fit_retrieval(pairs.TB(mask), ...
                    pairs.(operatorFields{operatorIndex})(mask));
                metricRows(end + 1, :) = {groupNames(groupIndex), ...
                    subsetNames(subsetIndex), operatorNames(operatorIndex), ...
                    fit.n, fit.correlation, fit.rSquared, fit.bias, fit.rmse, ...
                    fit.slope, fit.intercept, nnz(mask)}; %#ok<AGROW>
            end
            delta = pairs.TBSampleMinusParameter(mask);
            delta = delta(isfinite(delta));
            if isempty(delta)
                summary = nan(1, 8);
            else
                summary = [mean(delta), std(delta), median(delta), ...
                    prctile(delta, 5), prctile(delta, 95), ...
                    min(delta), max(delta), sqrt(mean(delta .^ 2))];
            end
            differenceRows(end + 1, :) = {groupNames(groupIndex), ...
                subsetNames(subsetIndex), numel(delta), summary(1), ...
                summary(2), summary(3), summary(4), summary(5), ...
                summary(6), summary(7), summary(8), nnz(mask)}; %#ok<AGROW>
        end
    end

    metrics = cell2table(metricRows, 'VariableNames', ...
        {'group','subset','operator','n','correlation','rSquared', ...
        'biasK','rmseK','slope','interceptK','selectedRows'});
    differences = cell2table(differenceRows, 'VariableNames', ...
        {'group','subset','n','meanK','stdK','medianK','p05K','p95K', ...
        'minK','maxK','rmsK','selectedRows'});
end
