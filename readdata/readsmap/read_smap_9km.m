function data = read_smap_9km(filePath, smapConfig)
%READ_SMAP_9KM 读取 SMAP SPL1CTB_E.004 9 km 北极半轨亮温。
%   DATA = READ_SMAP_9KM(FILEPATH) 返回文件元数据、全部前/后视观测
%   DATA.OBSERVATIONS，以及通过当前质量规则的 DATA.RECOMMENDED。
%
%   DATA = READ_SMAP_9KM(FILEPATH, SMAPCONFIG) 使用
%   RETRIEVAL_CONFIG().DATA.SMAP 中的亮温字段、质量位和有效范围设置。

    if nargin < 2 || isempty(smapConfig)
        config = retrieval_config();
        smapConfig = config.data.smap;
    end
    data = read_smap_l1c_product( ...
        filePath, smapConfig.enhanced, smapConfig);
end
