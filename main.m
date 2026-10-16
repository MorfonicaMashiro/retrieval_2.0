% 2017 年论文第 5.1 节：OIB 飞行日与同日 SMOS 配准。
% 在本目录运行 main。输出每个中心格的 TB、四个 OIB 均值和有效样本数。
% 反演须调用论文引用的完整 L 波段前向模型；此脚本不以近似模型伪造结果。

matlabRoot = fileparts(mfilename('fullpath'));
dataRoot = fullfile(fileparts(matlabRoot), 'data');
figureRoot = fullfile(matlabRoot, 'figures');
resultRoot = fullfile(matlabRoot, 'processed');
addpath(matlabRoot, fullfile(matlabRoot, 'readdata'), ...
    fullfile(matlabRoot, 'plot'));

gridVersion = "EASE"; % 原版 EASE 是显式假设；论文未给出足够参数区分版本。
useCache = true;
cacheFile = fullfile(resultRoot, 'pairs_2017_same_day_v2.mat');
if useCache && isfile(cacheFile)
    load(cacheFile, 'pairs', 'summary');
    fprintf('读取同日配对缓存：%d 个中心格。\n', height(pairs));
else
    timer = tic;
    pairs = table();
    summary = table();
    for year = [2012 2015]
        if year == 2012
            oibDir = fullfile(dataRoot, 'OIB', 'IDCSI4_L4', '2012');
            files = dir(fullfile(oibDir, 'IDCSI4_*.txt'));
        else
            oibDir = fullfile(dataRoot, 'OIB', ...
                'NSIDC-0708_QuickLook', '2015_GR_NASA');
            files = dir(fullfile(oibDir, 'OIB_*_IDCSI2.txt'));
        end
        [~, order] = sort({files.name});
        files = files(order);
        if isempty(files), error('main:NoOIB', '%d 年缺少 OIB 航次。', year); end
        for k = 1:numel(files)
            token = regexp(files(k).name, '(\d{8})', 'match', 'once');
            day = datetime(token, 'InputFormat', 'yyyyMMdd');
            oibFile = fullfile(files(k).folder, files(k).name);
            smosFile = fullfile(dataRoot, 'smos', num2str(year), ...
                ['SMOS_TB_north_' token '.nc']);
            oib = read_oib(oibFile);
            smos = read_smos(smosFile);
            smosEase = regrid_smos_to_ease(smos, gridVersion);
            [dailyPairs, stats] = collocate_2017_day(oib, smosEase, day);
            pairs = [pairs; dailyPairs]; %#ok<AGROW>
            summary = [summary; table(day, stats.rawPoints, ...
                stats.qualityPoints, stats.tbOverlapCells, ...
                stats.matchedCells, stats.mixedTypeCells, ...
                'VariableNames', {'day','rawPoints','qualityPoints', ...
                'tbOverlapCells','matchedCells','mixedTypeCells'})]; %#ok<AGROW>
            fprintf('%s：OIB 有效点 %d，匹配中心格 %d。\n', ...
                token, stats.qualityPoints, stats.matchedCells);
        end
    end
    if ~isfolder(resultRoot), mkdir(resultRoot); end
    save(cacheFile, 'pairs', 'summary', '-v7.3');
    writetable(pairs, fullfile(resultRoot, 'pairs_2017_same_day_v2.csv'));
    fprintf('全部航次同日配准：%.1f 秒。\n', toc(timer));
end

coverageFigure = fullfile(figureRoot, 'pairs_2017_coverage.png');
plot_2017_pairs(pairs, ease_grid(gridVersion), coverageFigure);
fprintf('共 %d 次飞行，%d 条中心格配对；位置图：%s\n', ...
    height(summary), height(pairs), coverageFigure);
% 工作区保留 pairs 与 summary。hi、hs 是 OIB 验证值，不是反演输出。
