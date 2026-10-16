function [pairs, stats] = collocate_2017_day(oib, smosEase, day)
%COLLOCATE_2017_DAY 按 2017 论文第 5.1 节汇集同日 OIB 与 SMOS。
%   每个有效 SMOS 中心格使用周围 3×3 EASE 格内的 OIB 原始记录。
%   对同一批有效记录分别求 FBsnow、hs、FBice、hi 的算术均值。
%   FBice = FBsnow - hs，是 OIB 几何冰面干舷代理。
%   本地质量字段中 0 表示无开阔水/薄冰/灰冰，my_ice_flag 为
%   0=FYI、1=MYI；九格内同时出现两种冰型的中心格会被排除。

    tbl = oib.table;
    required = {'pcnt_ow', 'pcnt_thin_ice', 'pcnt_grey_ice', 'my_ice_flag'};
    if ~all(ismember(required, tbl.Properties.VariableNames))
        error('collocate_2017_day:Fields', 'OIB 缺少冰型或混合像元字段。');
    end
    grid = smosEase.grid;
    fbSnow = double(tbl.mean_fb);
    hs = double(tbl.snow_depth);
    fbIce = fbSnow - hs;
    hi = double(tbl.thickness);
    iceType = double(tbl.my_ice_flag);

    % 论文要求排除混合冰型、开阔水和冰间水道；这里使用产品质量字段落实。
    % 论文第 5.1 节没有要求排除负冰面干舷，因此不额外设置该阈值。
    valid = isfinite(tbl.lat) & isfinite(tbl.lon) & ...
        isfinite(fbSnow) & isfinite(hs) & isfinite(hi) & ...
        isfinite(fbIce) & ...
        tbl.pcnt_ow == 0 & tbl.pcnt_thin_ice == 0 & ...
        tbl.pcnt_grey_ice == 0 & ismember(iceType, [0 1]);
    [row, col] = PSP2EASE(tbl.lat(valid), tbl.lon(valid), grid.version);
    located = isfinite(row) & isfinite(col);
    row = row(located);
    col = col(located);
    fbSnow = fbSnow(valid); fbSnow = fbSnow(located);
    hs = hs(valid); hs = hs(located);
    fbIce = fbIce(valid); fbIce = fbIce(located);
    hi = hi(valid); hi = hi(located);
    iceType = iceType(valid); iceType = iceType(located);

    shape = [grid.nRows, grid.nCols];
    nCells = prod(shape);
    index = sub2ind(shape, row, col);
    count = accumarray(index, 1, [nCells 1], @sum, 0);
    fyiCount = accumarray(index, double(iceType == 0), [nCells 1], @sum, 0);
    myiCount = accumarray(index, double(iceType == 1), [nCells 1], @sum, 0);
    fields = {fbSnow, hs, fbIce, hi};
    sums = cell(1, 4);
    kernel = ones(3);
    for k = 1:4
        cellSum = accumarray(index, fields{k}, [nCells 1], @sum, 0);
        sums{k} = conv2(reshape(cellSum, shape), kernel, 'same');
    end
    nOIB = conv2(reshape(count, shape), kernel, 'same');
    nFYI = conv2(reshape(fyiCount, shape), kernel, 'same');
    nMYI = conv2(reshape(myiCount, shape), kernel, 'same');
    good = nOIB > 0 & isfinite(smosEase.tb) & ...
        ((nFYI > 0 & nMYI == 0) | (nMYI > 0 & nFYI == 0));
    centers = find(good);
    [centerRow, centerCol] = ind2sub(shape, centers);
    sampleCount = nOIB(centers);
    [~, flightName] = fileparts(char(oib.source));
    n = numel(centers);
    pairs = table(repmat(day, n, 1), repmat(string(flightName), n, 1), ...
        centerRow, centerCol, smosEase.tb(centers), ...
        sums{1}(centers) ./ sampleCount, ...
        sums{2}(centers) ./ sampleCount, ...
        sums{3}(centers) ./ sampleCount, ...
        sums{4}(centers) ./ sampleCount, sampleCount, nMYI(centers) > 0, ...
        'VariableNames', {'day','flight','row','col','TB', ...
        'FBsnow','hs','FBice','hi','nOIB','isMYI'});
    stats = struct('rawPoints', height(tbl), 'qualityPoints', numel(index), ...
        'tbOverlapCells', nnz(nOIB > 0 & isfinite(smosEase.tb)), ...
        'matchedCells', n, 'mixedTypeCells', ...
        nnz(nOIB > 0 & nFYI > 0 & nMYI > 0 & isfinite(smosEase.tb)));
end
