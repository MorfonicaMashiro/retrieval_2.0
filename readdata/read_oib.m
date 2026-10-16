function data = read_oib(filePath)
%READ_OIB 读取 2012 年 IDCSI4 或 2015 年 NSIDC-0708 的 OIB 航次文本。
%   返回 source、table、snowFreeboard 和 iceFreeboardProxy。
%   后者为 mean_fb - snow_depth 推算量，不是独立雷达观测。

    filePath = char(string(filePath));
    if ~isfile(filePath)
        error('read_oib:MissingFile', 'OIB file does not exist: %s', filePath);
    end
    tbl = readtable(filePath, 'Delimiter', ',', ...
        'VariableNamingRule', 'preserve');
    required = {'lat', 'lon', 'mean_fb', 'snow_depth', 'thickness'};
    if ~all(ismember(required, tbl.Properties.VariableNames))
        error('read_oib:MissingVariable', ...
            'OIB file is missing one of: %s', strjoin(required, ', '));
    end

    % 两类本地 OIB 产品均用 -99999 表示数值字段缺测。
    for k = 1:width(tbl)
        values = tbl{:, k};
        if isnumeric(values)
            missing = values == -99999;
            if any(missing, 'all')
                values = double(values);
                values(missing) = NaN;
                tbl{:, k} = values;
            end
        end
    end

    data = struct();
    data.source = string(filePath);
    data.table = tbl;
    data.snowFreeboard = double(tbl.mean_fb);
    data.iceFreeboardProxy = data.snowFreeboard - double(tbl.snow_depth);
end
