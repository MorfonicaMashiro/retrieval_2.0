function [results, metricSummary, pairs, summary] = ...
        run_stride3_weighted_parameter_mean_inversion( ...
        targetDays, outputStem, maxPairsPerDay)
%RUN_STRIDE3_WEIGHTED_PARAMETER_MEAN_INVERSION 加权参数单柱联合反演。

    if nargin < 1, targetDays = datetime.empty(0,1); end
    if nargin < 3, maxPairsPerDay = Inf; end
    targetDays = dateshift(datetime(targetDays(:)), 'start', 'day');
    routeRoot = fileparts(mfilename('fullpath'));
    matlabRoot = fileparts(fileparts(fileparts(routeRoot)));
    dataRoot = fullfile(fileparts(matlabRoot), 'data');
    resultRoot = fullfile(routeRoot, 'results', 'inversion');
    figureRoot = fullfile(routeRoot, 'figures', 'inversion');
    addpath(matlabRoot, fullfile(matlabRoot,'readdata'), ...
        fullfile(matlabRoot,'plot'), routeRoot);
    if ~isfolder(resultRoot), mkdir(resultRoot); end
    if ~isfolder(figureRoot), mkdir(figureRoot); end

    retrievalConfig = retrieval_config();
    retrievalConfig.collocation.computeSampleMeanTB = false;
    retrievalConfig.collocation.parameterWeighting = ...
        "inverse_distance_to_center";
    retrievalConfig.collocation.distancePower = 1;
    model = retrievalConfig.model;
    operatorOptions = retrievalConfig.collocation;
    if nargin < 2 || strlength(string(outputStem)) == 0
        outputStem = "distance_weighted_parameter_inversion_" + ...
            retrievalConfig.output.variantTag + "_2012_2015";
    end
    outputStem = string(outputStem);
    flights = find_flight_files(dataRoot,targetDays);
    if isempty(flights)
        error('run_stride3_weighted_parameter_mean_inversion:NoFlights', ...
            '指定日期没有找到 OIB 航次。');
    end
    pairs = table(); summary = table();
    for flightIndex = 1:height(flights)
        day = flights.day(flightIndex);
        token = char(string(day,'yyyyMMdd'));
        oib = read_oib(char(flights.path(flightIndex)));
        smosPath = fullfile(dataRoot,'smos',char(string(year(day))), ...
            ['SMOS_TB_north_' token '.nc']);
        smosEase = grid_resampling.regrid_smos_to_ease( ...
            read_smos(smosPath), retrievalConfig.data.gridVersion);
        [dailyPairs,stats] = collocate_day( ...
            oib,smosEase,day,retrievalConfig);
        if isfinite(maxPairsPerDay) && height(dailyPairs)>maxPairsPerDay
            dailyPairs = dailyPairs(1:maxPairsPerDay,:);
        end
        pairs = [pairs;dailyPairs]; %#ok<AGROW>
        summary = [summary;table(day,height(dailyPairs),stats.rawPoints, ...
            stats.qualityPoints,'VariableNames', ...
            {'day','matchedPairs','rawPoints','qualityPoints'})]; %#ok<AGROW>
    end
    pairs = sortrows(pairs,{'day','row','col'});
    n = height(pairs);
    inputs = table(pairs.day,pairs.flight,pairs.row,pairs.col,pairs.TB, ...
        pairs.FBiceWeighted,pairs.FBsnowWeighted,pairs.isMYI, ...
        pairs.hiWeighted,pairs.hsWeighted,pairs.nOIB, ...
        pairs.nOIBWeightedSupport, ...
        'VariableNames',{'day','flight','row','col','TB','FBice', ...
        'FBsnow','isMYI','hi','hs','nOIB','nOIBWeightedSupport'});
    if retrievalConfig.temperature.source == "oib_kt19"
        inputs = addvars(inputs, pairs.surfaceTempCWeighted, ...
            'NewVariableNames', 'surfaceTempC');
    end
    hiRadar=nan(n,1); hsRadar=nan(n,1);
    hiLaser=nan(n,1); hsLaser=nan(n,1);
    hiLaserAlt=nan(n,1); hsLaserAlt=nan(n,1);
    nRadarRoots=zeros(n,1); nLaserRoots=zeros(n,1);
    tbFromOIB=pairs.TBOIBWeightedParameters;
    timer=tic;
    for pairIndex=1:n
        pairConfig = retrievalConfig;
        if retrievalConfig.temperature.source == "oib_kt19"
            pairConfig.model.surfaceTempC = inputs.surfaceTempC(pairIndex);
        end
        [radar,~]=invert_pair(inputs.TB(pairIndex), ...
            inputs.FBice(pairIndex),"ice",inputs.isMYI(pairIndex), ...
            pairConfig);
        [laser,~]=invert_pair(inputs.TB(pairIndex), ...
            inputs.FBsnow(pairIndex),"snow",inputs.isMYI(pairIndex), ...
            pairConfig);
        nRadarRoots(pairIndex)=numel(radar);
        nLaserRoots(pairIndex)=numel(laser);
        if ~isempty(radar)
            hiRadar(pairIndex)=radar(1).hi; hsRadar(pairIndex)=radar(1).hs;
        end
        if ~isempty(laser)
            hiLaser(pairIndex)=laser(1).hi; hsLaser(pairIndex)=laser(1).hs;
        end
        if numel(laser)>1
            hiLaserAlt(pairIndex)=laser(2).hi;
            hsLaserAlt(pairIndex)=laser(2).hs;
        end
    end
    results=addvars(inputs,hiRadar,hsRadar,hiLaser,hsLaser, ...
        hiLaserAlt,hsLaserAlt,nRadarRoots,nLaserRoots,tbFromOIB);
    nearTB=abs(results.tbFromOIB-results.TB)<=1.5;
    sortedCounts=sort(results.nOIB(isfinite(results.nOIB)));
    coverageThreshold=sortedCounts( ...
        max(1,ceil(0.95*numel(sortedCounts))));
    highCoverage=results.nOIB>=coverageThreshold;
    metricSummary=summarize_retrieval_subset( ...
        results,true(n,1),nearTB);
    highCoverageFit=fit_retrieval( ...
        results.TB(highCoverage),results.tbFromOIB(highCoverage));
    highCoverageRow=metricSummary(1,:);
    highCoverageRow.metric="TB：OIB 样本数前 5%";
    fitFields={'n','correlation','rSquared','bias','rmse','slope','intercept'};
    for fieldIndex=1:numel(fitFields)
        fieldName=fitFields{fieldIndex};
        highCoverageRow.(fieldName)=highCoverageFit.(fieldName);
    end
    metricSummary=[metricSummary(1,:);highCoverageRow; ...
        metricSummary(2:end,:)];
    metricByYear=summarize_retrieval_by_year(results,nearTB, ...
        fullfile(resultRoot,outputStem+"_metrics_by_year.csv"));
    outputFile=fullfile(resultRoot,outputStem+".mat");
    save(outputFile,'pairs','summary','results','metricSummary', ...
        'metricByYear','nearTB','highCoverage','coverageThreshold', ...
        'model','operatorOptions','retrievalConfig','-v7.3');
    writetable(results,fullfile(resultRoot,outputStem+"_retrieval.csv"));
    writetable(summary,fullfile(resultRoot,outputStem+"_summary.csv"));
    writetable(metricSummary,fullfile(resultRoot,outputStem+"_metrics.csv"));
    plot_retrieval_validation(results,nearTB,highCoverage,figureRoot, ...
        '到 SMOS 中心反距离加权参数→单一冰雪柱→单次 TB');
    fprintf('完成 %d 条加权参数反演，用时 %.1f min：%s\n', ...
        n,toc(timer)/60,outputFile);
end

function flights=find_flight_files(dataRoot,targetDays)
    specs={fullfile(dataRoot,'OIB','IDCSI4_L4','2012'),'IDCSI4_*.txt'; ...
        fullfile(dataRoot,'OIB','NSIDC-0708_QuickLook','2015_GR_NASA'), ...
        'OIB_*_IDCSI2*.txt'};
    days=datetime.empty(0,1); paths=strings(0,1);
    for specIndex=1:size(specs,1)
        files=dir(fullfile(specs{specIndex,1},specs{specIndex,2}));
        for fileIndex=1:numel(files)
            token=regexp(files(fileIndex).name,'(\d{8})','match','once');
            day=datetime(token,'InputFormat','yyyyMMdd');
            if isempty(targetDays)||ismember(day,targetDays)
                days(end+1,1)=day; %#ok<AGROW>
                paths(end+1,1)=string(fullfile( ...
                    files(fileIndex).folder,files(fileIndex).name)); %#ok<AGROW>
            end
        end
    end
    flights=sortrows(table(days,paths, ...
        'VariableNames',{'day','path'}),'day');
end
