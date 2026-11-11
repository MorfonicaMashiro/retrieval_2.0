function [tailRows, qualitySummary] = inspect_stride3_tail_smos_quality()
%INSPECT_STRIDE3_TAIL_SMOS_QUALITY 检查长尾点所在 SMOS 源像元质量字段。
%   只读取 |TB_sim-TB_SMOS|>=10 K 所在的飞行日，并将同日其他配对作为
%   对照，输出 nPair、TB_uncertainty 和 RFI_ratio。

    experimentRoot = fileparts(mfilename('fullpath'));
    matlabRoot = fileparts(fileparts(fileparts(experimentRoot)));
    workspaceRoot = fileparts(matlabRoot);
    resultRoot = fullfile(experimentRoot, 'results');
    dataRoot = fullfile(workspaceRoot, 'data');
    addpath(matlabRoot, fullfile(matlabRoot, 'readdata'));

    sourceFile = fullfile(resultRoot, ...
        'stride3_parameter_mean_2012_2015_retrieval.csv');
    results = readtable(sourceFile);
    results.differenceK = results.tbFromOIB - results.TB;
    results.absDifferenceK = abs(results.differenceK);
    tailThresholdK = 10;
    tailDays = unique(results.day(results.absDifferenceK >= tailThresholdK));
    selected = ismember(results.day, tailDays);
    inspected = results(selected,:);
    inspected.isTail10K = inspected.absDifferenceK >= tailThresholdK;
    inspected.smosNPair = nan(height(inspected),1);
    inspected.smosTBUncertaintyK = nan(height(inspected),1);
    inspected.smosRFIRatio = nan(height(inspected),1);
    inspected.regriddedTBCheckK = nan(height(inspected),1);
    oldEaseGrid = grid_resampling.ease_grid("EASE");
    centerX = oldEaseGrid.xCenters(inspected.col).';
    centerY = oldEaseGrid.yCenters(inspected.row).';
    [inspected.centerLatitude, inspected.centerLongitude] = projinv( ...
        oldEaseGrid.projection, centerX, centerY);

    for dayIndex = 1:numel(tailDays)
        currentDay = tailDays(dayIndex);
        token = char(string(currentDay, 'yyyyMMdd'));
        smosFile = fullfile(dataRoot, 'smos', ...
            char(string(year(currentDay))), ['SMOS_TB_north_' token '.nc']);
        smosEase = grid_resampling.regrid_smos_to_ease( ...
            read_smos(smosFile), "EASE");
        rows = find(inspected.day == currentDay);
        linear = sub2ind(size(smosEase.tb), ...
            inspected.row(rows), inspected.col(rows));
        inspected.smosNPair(rows) = smosEase.nPair(linear);
        inspected.smosTBUncertaintyK(rows) = ...
            smosEase.tbUncertainty(linear);
        inspected.smosRFIRatio(rows) = smosEase.rfiRatio(linear);
        inspected.regriddedTBCheckK(rows) = smosEase.tb(linear);
    end
    maximumTBMismatchK = max(abs( ...
        inspected.regriddedTBCheckK - inspected.TB), [], 'omitnan');
    if maximumTBMismatchK > 1e-9
        error('inspect_stride3_tail_smos_quality:TBMismatch', ...
            '重新读取的 SMOS TB 与冻结结果不一致，最大差值 %.3g K。', ...
            maximumTBMismatchK);
    end

    tailRows = inspected(inspected.isTail10K,:);
    tailRows = sortrows(tailRows, 'absDifferenceK', 'descend');
    qualitySummary = [quality_row("同日非长尾", ...
        inspected(~inspected.isTail10K,:)); ...
        quality_row("|差值|≥10 K", inspected(inspected.isTail10K,:))];

    writetable(tailRows, fullfile(resultRoot, ...
        'stride3_tb_residual_abs_ge_10K_with_smos_quality.csv'));
    writetable(qualitySummary, fullfile(resultRoot, ...
        'stride3_tb_tail_smos_quality_summary.csv'));
    fprintf('重新读取 TB 最大差值：%.3g K\n', maximumTBMismatchK);
    disp(qualitySummary);
end


function output = quality_row(label, rows)
    output = table(string(label), height(rows), ...
        mean(rows.differenceK, 'omitnan'), ...
        median(rows.smosNPair, 'omitnan'), ...
        median(rows.smosTBUncertaintyK, 'omitnan'), ...
        prctile(rows.smosTBUncertaintyK, 95), ...
        median(rows.smosRFIRatio, 'omitnan'), ...
        prctile(rows.smosRFIRatio, 95), ...
        'VariableNames', {'group','n','meanDifferenceK', ...
        'medianSMOSNPair','medianTBUncertaintyK', ...
        'percentile95TBUncertaintyK','medianRFIRatio', ...
        'percentile95RFIRatio'});
end
