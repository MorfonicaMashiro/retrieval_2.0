% 2017 年论文第 5.1 节：OIB 飞行日与同日 SMOS 配准。
% 在本目录运行 main。输出每个中心格的 TB、四个 OIB 均值和有效样本数。
% 反演须调用论文引用的完整 L 波段前向模型；此脚本不以近似模型伪造结果。
clc
clear all


matlabRoot = fileparts(mfilename('fullpath'));
dataRoot = fullfile(fileparts(matlabRoot), 'data');
figureRoot = fullfile(matlabRoot, 'figures');
resultRoot = fullfile(matlabRoot, 'processed');
addpath(matlabRoot, fullfile(matlabRoot, 'readdata'), ...
    fullfile(matlabRoot, 'plot'));

gridVersion = "EASE"; % 原版 EASE 是显式假设；论文未给出足够参数区分版本。
yearsToProcess = [2012,2015];
useCache = true;
yearTag = strjoin(string(yearsToProcess), '_');
cacheStem = "pairs_" + yearTag + "_same_day";
cacheFile = fullfile(resultRoot, cacheStem + ".mat");
cacheIsValid = false;
if useCache && isfile(cacheFile)
    cacheData = load(cacheFile, 'pairs', 'summary');
    if isfield(cacheData, 'pairs') && isfield(cacheData, 'summary')
        cachedYears = unique(year(cacheData.summary.day)).';
        cacheIsValid = isequal(cachedYears, yearsToProcess);
    end
    if ~cacheIsValid
        warning('main:IncompleteCache', ...
            '缓存 %s 年份不符或内容不完整，将按 %s 重建。', ...
            cacheFile, strjoin(string(yearsToProcess), ','));
    end
end
if cacheIsValid
    pairs = cacheData.pairs;
    summary = cacheData.summary;
    fprintf('读取 %s 同日配对缓存：%d 个中心格。\n', ...
        strjoin(string(yearsToProcess), ','), height(pairs));
else
    timer = tic;
    pairs = table();
    summary = table();
    flightFiles = cell(size(yearsToProcess));
    totalFlights = 0;
    for yearIndex = 1:numel(yearsToProcess)
        currentYear = yearsToProcess(yearIndex);
        if currentYear <= 2013
            % 2012–2013 必须使用正式 IDCSI4，不回退到 Quick Look。
            oibDir = fullfile(dataRoot, 'OIB', 'IDCSI4_L4', ...
                num2str(currentYear));
            files = dir(fullfile(oibDir, 'IDCSI4_*.txt'));
        else
            oibDir = fullfile(dataRoot, 'OIB', ...
                'NSIDC-0708_QuickLook', sprintf('%d_GR_NASA', currentYear));
            files = dir(fullfile(oibDir, 'OIB_*_IDCSI2*.txt'));
        end
        [~, order] = sort({files.name});
        files = files(order);
        if isempty(files)
            error('main:NoOIB', ...
                '%d 年缺少指定 OIB 产品：%s。', currentYear, oibDir);
        end
        flightFiles{yearIndex} = files;
        totalFlights = totalFlights + numel(files);
    end
    completedFlights = 0;
    for yearIndex = 1:numel(yearsToProcess)
        currentYear = yearsToProcess(yearIndex);
        files = flightFiles{yearIndex};
        for k = 1:numel(files)
            token = regexp(files(k).name, '(\d{8})', 'match', 'once');
            day = datetime(token, 'InputFormat', 'yyyyMMdd');
            oibFile = fullfile(files(k).folder, files(k).name);
            smosFile = fullfile(dataRoot, 'smos', num2str(currentYear), ...
                ['SMOS_TB_north_' token '.nc']);
            if ~isfile(smosFile)
                error('main:MissingSMOS', ...
                    '%s 缺少同日 SMOS 文件：%s', token, smosFile);
            end
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
            completedFlights = completedFlights + 1;
            elapsedSeconds = toc(timer);
            etaSeconds = elapsedSeconds / completedFlights * ...
                (totalFlights - completedFlights);
            fprintf(['[%d/%d, %5.1f%%] %s：OIB 有效点 %d，' ...
                '匹配中心格 %d；已用 %.1f min，预计剩余 %.1f min。\n'], ...
                completedFlights, totalFlights, ...
                100 * completedFlights / totalFlights, token, ...
                stats.qualityPoints, stats.matchedCells, ...
                elapsedSeconds / 60, etaSeconds / 60);
        end
    end
    if ~isfolder(resultRoot), mkdir(resultRoot); end
    save(cacheFile, 'pairs', 'summary', '-v7.3');
    writetable(pairs, fullfile(resultRoot, cacheStem + ".csv"));
    fprintf('全部航次同日配准：%.1f 秒。\n', toc(timer));
end

coverageFigure = fullfile(figureRoot, ...
    "pairs_" + yearTag + "_coverage.png");
plot_2017_pairs(pairs, ease_grid(gridVersion), coverageFigure);
fprintf('共 %d 次飞行，%d 条中心格配对；位置图：%s\n', ...
    height(summary), height(pairs), coverageFigure);
% 工作区保留 pairs 与 summary。hi、hs 是 OIB 验证值，不是反演输出。
