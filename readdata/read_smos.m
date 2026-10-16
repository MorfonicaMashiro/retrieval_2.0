function data = read_smos(filePath)
%READ_SMOS 读取一份 ICDC SMOS L3B 北半球逐日 NetCDF 文件。
%   输出矩阵按 [行,列] 排列，大小为 [896,608]，便于直接绘图。
%   MATLAB ncread 按 [x,y,time] 读出；去掉单元素 time 维后需转置。
%   保留原生 NSIDC 极地立体投影网格，不做 EASE 转换。

    filePath = char(string(filePath));
    if ~isfile(filePath)
        error('read_smos:MissingFile', 'SMOS file does not exist: %s', filePath);
    end
    [~, stem] = fileparts(filePath);
    token = regexp(stem, '^SMOS_TB_north_(\d{8})$', 'tokens', 'once');
    if isempty(token)
        error('read_smos:Filename', 'Unexpected SMOS filename: %s', stem);
    end

    data = struct();
    data.source = string(filePath);
    data.day = datetime(token{1}, 'InputFormat', 'yyyyMMdd');
    data.tb = read_field(filePath, 'TB');
    data.latitude = read_field(filePath, 'latitude');
    data.longitude = read_field(filePath, 'longitude');
    data.tbUncertainty = read_field(filePath, 'TB_uncertainty');
    data.nPair = read_field(filePath, 'nPair');
    data.rfiRatio = read_field(filePath, 'RFI_ratio');

    fields = {'latitude', 'longitude', 'tbUncertainty', 'nPair', 'rfiRatio'};
    for k = 1:numel(fields)
        if ~isequal(size(data.tb), size(data.(fields{k})))
            error('read_smos:SizeMismatch', ...
                'SMOS field %s has a different raster size.', fields{k});
        end
    end
end

function values = read_field(filePath, variableName)
    values = double(squeeze(ncread(filePath, variableName))).';
    % 此批文件的数值缺测标记为 -999，读取后统一转为 NaN。
    values(values == -999) = NaN;
end
