function tests = test_read_smap
%TEST_READ_SMAP SMAP 36 km/9 km 正式读取接口的结构回归测试。
    tests = functiontests(localfunctions);
end

function testReadStandard36Km(testCase)
    fixture = create_fixture(testCase, "36km");
    data = read_smap_36km(fixture.filePath);

    verifyEqual(testCase, data.product, "SPL1CTB.006");
    verifyEqual(testCase, data.resolutionKm, 36);
    verifyEqual(testCase, data.northGroup, "/North_Polar_Projection");
    verifyEqual(testCase, height(data.observations), 6);
    verifyEqual(testCase, height(data.recommended), 2);
    verifyEqual(testCase, unique(data.observations.look), ["aft"; "fore"]);
    verifyTrue(testCase, all(isnat(data.observations.timeUtc)));
end

function testReadEnhanced9Km(testCase)
    fixture = create_fixture(testCase, "9km");
    data = read_smap_9km(fixture.filePath);

    verifyEqual(testCase, data.product, "SPL1CTB_E.004");
    verifyEqual(testCase, data.resolutionKm, 9);
    verifyEqual(testCase, data.northGroup, ...
        "/EASE2_North_Polar_Projection");
    verifyEqual(testCase, height(data.observations), 6);
    verifyEqual(testCase, height(data.recommended), 2);
    verifyEqual(testCase, data.orbitNumber, 1);
    verifyEqual(testCase, data.orbitDirection, "A");
end

function testSurfaceCorrectedRequestRequiresField(testCase)
    fixture = create_fixture(testCase, "36km");
    config = retrieval_config();
    config.data.smap.tbField = "surface_corrected";
    verifyError(testCase, ...
        @() read_smap_36km(fixture.filePath, config.data.smap), ...
        'read_smap_l1c_product:MissingSurfaceCorrectedTB');
end

function fixture = create_fixture(testCase, product)
    folder = tempname;
    mkdir(folder);
    testCase.addTeardown(@() remove_fixture(folder));
    if product == "36km"
        filename = 'SMAP_L1C_TB_00001_A_20160419T000000_R19240_001.h5';
        group = "/North_Polar_Projection";
    else
        filename = 'SMAP_L1C_TB_E_00001_A_20160419T000000_R19240_001.h5';
        group = "/EASE2_North_Polar_Projection";
    end
    filePath = fullfile(folder, filename);

    write_dataset(filePath, group + "/cell_lat", [70; 71; 72], 'double');
    write_dataset(filePath, group + "/cell_lon", [-150; -149; -148], 'double');
    write_dataset(filePath, group + "/cell_row", uint16([1; 2; 3]), 'uint16');
    write_dataset(filePath, group + "/cell_column", ...
        uint16([4; 5; 6]), 'uint16');
    for look = ["fore", "aft"]
        suffix = "_" + look;
        write_dataset(filePath, group + "/cell_tb_h" + suffix, ...
            single([70; 250; 320]), 'single');
        write_dataset(filePath, group + "/cell_tb_v" + suffix, ...
            single([80; 260; 300]), 'single');
        write_dataset(filePath, group + "/cell_tb_qual_flag_h" + suffix, ...
            uint16([0; 1; 0]), 'uint16');
        write_dataset(filePath, group + "/cell_tb_qual_flag_v" + suffix, ...
            uint16([0; 0; 0]), 'uint16');
    end
    fixture = struct('filePath', filePath);
end

function write_dataset(filePath, datasetPath, values, dataType)
    h5create(filePath, char(datasetPath), size(values), 'Datatype', dataType);
    h5write(filePath, char(datasetPath), values);
end

function remove_fixture(folder)
    if isfolder(folder)
        rmdir(folder, 's');
    end
end
