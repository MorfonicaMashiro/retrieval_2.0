function metricByYear = summarize_retrieval_by_year( ...
        results, nearTB, outputCsv, nearTbThresholdK)
%SUMMARIZE_RETRIEVAL_BY_YEAR 逐年汇总亮温和两条联合反演路线的验证指标。
%   “all”使用该年全部有效解；“tb_within_1p5K”进一步要求由 OIB
%   冰厚/雪深模拟的亮温与 SMOS 观测亮温之差不超过 1.5 K。

    arguments
        results table
        nearTB (:,1) logical
        outputCsv {mustBeTextScalar} = ""
        nearTbThresholdK (1,1) double {mustBePositive} = 1.5
    end
    if height(results) ~= numel(nearTB)
        error('summarize_retrieval_by_year:MaskSize', ...
            'nearTB 与结果表行数不一致。');
    end

    years = unique(year(results.day));
    nearPrefix = "|TB误差|<=" + string(nearTbThresholdK) + " K：";
    metricNames = ["TB：全部配对"; ...
        "TB+FBice：冰厚"; "TB+FBice：雪深"; ...
        "TB+FBsnow：冰厚"; "TB+FBsnow：雪深"; ...
        nearPrefix + "TB+FBice 冰厚"; ...
        nearPrefix + "TB+FBice 雪深"; ...
        nearPrefix + "TB+FBsnow 冰厚"; ...
        nearPrefix + "TB+FBsnow 雪深"];
    nMetrics = numel(metricNames);
    nRows = numel(years) * nMetrics;
    yearColumn = zeros(nRows, 1);
    metricColumn = strings(nRows, 1);
    nColumn = zeros(nRows, 1);
    correlation = nan(nRows, 1);
    rSquared = nan(nRows, 1);
    bias = nan(nRows, 1);
    rmse = nan(nRows, 1);
    slope = nan(nRows, 1);
    intercept = nan(nRows, 1);

    row = 0;
    for yearValue = reshape(years, 1, [])
        yearMask = year(results.day) == yearValue;
        restricted = yearMask & nearTB;
        fits = { ...
            fit_retrieval(results.TB(yearMask), results.tbFromOIB(yearMask)); ...
            fit_retrieval(results.hi(yearMask), results.hiRadar(yearMask)); ...
            fit_retrieval(results.hs(yearMask), results.hsRadar(yearMask)); ...
            fit_retrieval(results.hi(yearMask), results.hiLaser(yearMask)); ...
            fit_retrieval(results.hs(yearMask), results.hsLaser(yearMask)); ...
            fit_retrieval(results.hi(restricted), results.hiRadar(restricted)); ...
            fit_retrieval(results.hs(restricted), results.hsRadar(restricted)); ...
            fit_retrieval(results.hi(restricted), results.hiLaser(restricted)); ...
            fit_retrieval(results.hs(restricted), results.hsLaser(restricted))};
        for metricIndex = 1:nMetrics
            row = row + 1;
            currentFit = fits{metricIndex};
            yearColumn(row) = yearValue;
            metricColumn(row) = metricNames(metricIndex);
            nColumn(row) = currentFit.n;
            correlation(row) = currentFit.correlation;
            rSquared(row) = currentFit.rSquared;
            bias(row) = currentFit.bias;
            rmse(row) = currentFit.rmse;
            slope(row) = currentFit.slope;
            intercept(row) = currentFit.intercept;
        end
    end
    metricByYear = table(yearColumn, metricColumn, nColumn, correlation, ...
        rSquared, bias, rmse, slope, intercept, 'VariableNames', ...
        {'year','metric','n','correlation','rSquared','bias','rmse', ...
        'slope','intercept'});

    outputCsv = char(string(outputCsv));
    if ~isempty(outputCsv)
        writetable(metricByYear, outputCsv);
    end
end
