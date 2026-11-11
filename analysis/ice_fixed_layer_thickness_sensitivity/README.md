# 固定冰层厚度敏感性测试

本目录是独立测试分支，不修改冻结的 `retrieval_config.m`，也不改变
`tb_forward.m` 默认的 20 层冰、1 层雪。

测试论文表 1 的八个理想场景，冰层目标厚度分别为 1、5、10 cm。由于这些
场景的总冰厚都可被三种目标厚度严格整除，使用

\[
N_{ice}=h_i/\Delta z
\]

即可调用原前向模型得到严格等厚的冰层，不需要在最后设置残余薄层。物理口径
保持 Burke 递推、表温 −30 °C、雪 1 层、不启用 Vb 70‰ 截断。

运行：

```matlab
cd('D:\海冰厚度联合反演\retrieval_matlab_ver')
addpath('analysis\ice_fixed_layer_thickness_sensitivity')
run_fixed_ice_layer_thickness_sensitivity
```

结果写入本目录下的 `results/`，不会覆盖基线验证或反演结果。

## 为什么暂不细分雪层

当前干雪模型采用统一密度、稳态线性温度廓线和 Tiuri 介电关系。同类雪层之间
在 Burke 递推中不设置内部 Fresnel 反射，因此增加雪层数主要只是更细地积分
温度引起的微弱介电损耗变化。L 波段下干雪本身损耗很小，论文复现配置也明确
采用 1 层雪。若没有逐层密度、湿度、冰壳或粒径资料，增加雪层数并不会引入
新的雪物理，只会增加计算量。
