# MATLAB 测试文件

这里集中存放项目的自动测试代码。日常运行和修改主程序时，不需要打开或改写这些文件。

在 MATLAB 中进入项目目录后，只运行：

```matlab
cd('D:/海冰厚度联合反演/retrieval_matlab_ver');
report = run_regression_suite();
```

根目录的 `run_regression_suite.m` 只是一个简短入口，实际测试均在本文件夹中：

| 文件 | 检查内容 |
| --- | --- |
| `test_retrieval_config.m` | 默认配置、覆盖项、非法配置和自动输出命名 |
| `test_tb_forward.m` | 前向模型公式、温度、Vb、分层、数值范围和论文场景 |
| `test_inversion_closure.m` | 8 个理想场景、两条干舷路线的反演自闭环 |
| `test_read_smap.m` | SMAP 36/9 km HDF5 字段、质量掩膜和接口 |
| `test_regrid_to_ease2_12p5.m` | EASE2-N 12.5 km 网格定义和记录定位 |
| `test_smap36_ease_overlap_weights.m` | SMAP 36 km 与原版 EASE 像元的面积闭合 |
| `run_regression_suite_impl.m` | 汇总上述测试、4 个算子验证及 MATLAB 静态检查 |

这些测试使用固定且可核对的输入，不是随机输入。它们主要保证代码行为没有被意外改坏，但不能代替 2012/2015 真实数据的全量重算和科学结果比较。
