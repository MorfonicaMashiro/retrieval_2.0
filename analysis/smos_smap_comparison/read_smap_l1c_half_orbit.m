function observations = read_smap_l1c_half_orbit(filePath)
%READ_SMAP_L1C_HALF_ORBIT 读取 SPL1CTB V006 北极投影半轨亮温。
%   只保留主质量位(bit 0)为 0 且 H/V 同时有效的前、后视观测。
%   输出为长表，每行是一条网格观测，不在此处做日平均。

    filePath = char(string(filePath));
    if ~isfile(filePath)
        error('read_smap_l1c_half_orbit:MissingFile', ...
            'SMAP file does not exist: %s', filePath);
    end
    group = locate_north_group(filePath);
    latitude = read_vector(filePath, group + "/cell_lat");
    longitude = read_vector(filePath, group + "/cell_lon");
    if numel(latitude) ~= numel(longitude)
        error('read_smap_l1c_half_orbit:GeolocationSize', ...
            'cell_lat and cell_lon have different sizes: %s', filePath);
    end

    observations = table();
    looks = ["fore", "aft"];
    for lookIndex = 1:numel(looks)
        look = looks(lookIndex);
        tbH = read_vector(filePath, group + "/cell_tb_h_" + look);
        tbV = read_vector(filePath, group + "/cell_tb_v_" + look);
        flagH = read_vector(filePath, ...
            group + "/cell_tb_qual_flag_h_" + look);
        flagV = read_vector(filePath, ...
            group + "/cell_tb_qual_flag_v_" + look);
        if any([numel(tbH), numel(tbV), numel(flagH), numel(flagV)] ...
                ~= numel(latitude))
            error('read_smap_l1c_half_orbit:FieldSize', ...
                'SMAP fields have inconsistent sizes: %s', filePath);
        end
        acceptableH = bitget(uint16(flagH), 1) == 0;
        acceptableV = bitget(uint16(flagV), 1) == 0;
        valid = isfinite(latitude) & isfinite(longitude) & ...
            latitude >= 0 & abs(longitude) <= 180 & ...
            isfinite(tbH) & isfinite(tbV) & ...
            tbH >= 50 & tbH <= 350 & tbV >= 50 & tbV <= 350 & ...
            acceptableH & acceptableV;
        n = nnz(valid);
        current = table(latitude(valid), longitude(valid), ...
            tbH(valid), tbV(valid), repmat(look, n, 1), ...
            'VariableNames', {'latitude','longitude','tbH','tbV','look'});
        observations = [observations; current]; %#ok<AGROW>
    end
    [~, name, extension] = fileparts(filePath);
    observations.sourceFile = repmat(string(name) + string(extension), ...
        height(observations), 1);
end

function group = locate_north_group(filePath)
    info = h5info(filePath);
    groupNames = string({info.Groups.Name});
    candidates = ["/North_Polar_Projection", ...
        "/EASE2_North_Polar_Projection"];
    selected = candidates(ismember(candidates, groupNames));
    if isempty(selected)
        error('read_smap_l1c_half_orbit:NorthGroup', ...
            'Could not find the north polar projection group in %s. Groups: %s', ...
            filePath, strjoin(groupNames, ', '));
    end
    group = selected(1);
end

function values = read_vector(filePath, datasetPath)
    values = double(h5read(filePath, char(datasetPath)));
    values = values(:);
    values(values <= -9990) = NaN;
end
