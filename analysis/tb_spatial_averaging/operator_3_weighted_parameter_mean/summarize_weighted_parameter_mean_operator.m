function [metrics, differences] = ...
        summarize_weighted_parameter_mean_operator(pairs)
%SUMMARIZE_WEIGHTED_PARAMETER_MEAN_OPERATOR 汇总算术与加权参数算子。

    required = {'day','TB','nOIB','nOIBWeightedSupport', ...
        'TBOIBMeanParameters','TBOIBWeightedParameters', ...
        'TBWeightedMinusArithmetic','hi','hs','hiWeighted','hsWeighted'};
    if ~all(ismember(required, pairs.Properties.VariableNames))
        error('summarize_weighted_parameter_mean_operator:Fields', ...
            '输入配对表缺少加权参数算子字段。');
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
        groupNames(end + 1, 1) = ...
            "day_" + string(currentDay, 'yyyyMMdd'); %#ok<AGROW>
        groupMasks{end + 1, 1} = pairs.day == currentDay; %#ok<AGROW>
    end

    subsetNames = ["all", "nOIB_ge_50"];
    operatorNames = ["mean_parameters", "weighted_parameters"];
    operatorFields = {'TBOIBMeanParameters', 'TBOIBWeightedParameters'};
    metricRows = cell(0, 11);
    differenceRows = cell(0, 13);
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
            deltaTB = pairs.TBWeightedMinusArithmetic(mask);
            deltaHi = pairs.hiWeighted(mask) - pairs.hi(mask);
            deltaHs = pairs.hsWeighted(mask) - pairs.hs(mask);
            finite = isfinite(deltaTB) & isfinite(deltaHi) & isfinite(deltaHs);
            deltaTB = deltaTB(finite);
            deltaHi = deltaHi(finite);
            deltaHs = deltaHs(finite);
            if isempty(deltaTB)
                summary = nan(1, 8);
                meanDeltaHi = NaN;
                meanDeltaHs = NaN;
            else
                summary = [mean(deltaTB), std(deltaTB), median(deltaTB), ...
                    prctile(deltaTB, 5), prctile(deltaTB, 95), ...
                    min(deltaTB), max(deltaTB), sqrt(mean(deltaTB .^ 2))];
                meanDeltaHi = mean(deltaHi);
                meanDeltaHs = mean(deltaHs);
            end
            differenceRows(end + 1, :) = {groupNames(groupIndex), ...
                subsetNames(subsetIndex), numel(deltaTB), summary(1), ...
                summary(2), summary(3), summary(4), summary(5), ...
                summary(6), summary(7), summary(8), meanDeltaHi, ...
                meanDeltaHs}; %#ok<AGROW>
        end
    end

    metrics = cell2table(metricRows, 'VariableNames', ...
        {'group','subset','operator','n','correlation','rSquared', ...
        'biasK','rmseK','slope','interceptK','selectedRows'});
    differences = cell2table(differenceRows, 'VariableNames', ...
        {'group','subset','n','meanTBK','stdTBK','medianTBK', ...
        'p05TBK','p95TBK','minTBK','maxTBK','rmsTBK', ...
        'meanDeltaHiM','meanDeltaHsM'});
end
