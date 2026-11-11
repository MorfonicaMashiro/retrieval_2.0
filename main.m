% OIB 飞行日与同日网格化 SMOS 配准入口。
% 在本目录运行 main。输出每个中心格的 TB、四个 OIB 均值和有效样本数。
% 反演须调用论文引用的完整 L 波段前向模型；此脚本不以近似模型伪造结果。
clc
clearvars -except useCache retrievalConfig


matlabRoot = fileparts(mfilename('fullpath'));
addpath(matlabRoot);
if ~exist('retrievalConfig', 'var') || isempty(retrievalConfig)
    retrievalConfig = retrieval_config();
end
addpath(retrievalConfig.paths.readDataRoot, ...
    retrievalConfig.paths.readSmapRoot, retrievalConfig.paths.plotRoot);
dataRoot = retrievalConfig.paths.dataRoot;
figureRoot = retrievalConfig.paths.figureRoot;
resultRoot = retrievalConfig.paths.processedRoot;
gridVersion = retrievalConfig.data.gridVersion;
yearsToProcess = retrievalConfig.data.years;
if retrievalConfig.data.brightnessTemperatureSource ~= "smos"
    error('main:SMAPPipelineNotImplemented', ...
        ['SMAP 读取接口和配置已就绪，但半轨逐格 UTC 选择与 OIB 配准尚未接入 ', ...
        'main.m。请勿在接入前把 SMAP 配置误当成已运行的主结果。']);
end
if ~exist('useCache', 'var'), useCache = retrievalConfig.runtime.useCache; end
cacheStem = retrievalConfig.output.pairsStem;
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
        oibConfig = retrievalConfig.data.oib;
        if currentYear <= oibConfig.formalLastYear
            % 2012–2013 必须使用正式 IDCSI4，不回退到 Quick Look。
            oibDir = fullfile(dataRoot, ...
                oibConfig.formalRelativeDirectory, ...
                num2str(currentYear));
            files = dir(fullfile(oibDir, oibConfig.formalFilePattern));
        else
            quickLookFolder = sprintf('%d%s', currentYear, ...
                oibConfig.quickLookYearSuffix);
            oibDir = fullfile(dataRoot, ...
                oibConfig.quickLookRelativeDirectory, quickLookFolder);
            files = dir(fullfile(oibDir, oibConfig.quickLookFilePattern));
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
            smosConfig = retrievalConfig.data.smos;
            smosFile = fullfile(dataRoot, smosConfig.relativeDirectory, ...
                num2str(currentYear), sprintf(smosConfig.filePattern, token));
            if ~isfile(smosFile)
                error('main:MissingSMOS', ...
                    '%s 缺少同日 SMOS 文件：%s', token, smosFile);
            end
            oib = read_oib(oibFile);
            smos = read_smos(smosFile);
            smosEase = grid_resampling.regrid_smos_to_ease( ...
                smos, gridVersion);
            [dailyPairs, stats] = collocate_day( ...
                oib, smosEase, day, retrievalConfig);
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
    retrievalConfig.output.coverageFigureStem + ".png");
plot_collocated_pairs(pairs, ...
    grid_resampling.ease_grid(gridVersion), coverageFigure);
fprintf('共 %d 次飞行，%d 条中心格配对；位置图：%s\n', ...
    height(summary), height(pairs), coverageFigure);
% 工作区保留 pairs 与 summary。hi、hs 是 OIB 验证值，不是反演输出。
