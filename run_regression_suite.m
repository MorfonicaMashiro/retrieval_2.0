%回归测试运行命令
function report = run_regression_suite()
%RUN_REGRESSION_SUITE 一键运行 tests 文件夹中的完整回归测试。

    matlabRoot = fileparts(mfilename('fullpath'));
    addpath(fullfile(matlabRoot, 'tests'));
    report = run_regression_suite_impl();
end
