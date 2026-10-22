function fit = fit_retrieval_2017(observed, retrieved)
%FIT_RETRIEVAL_2017 计算论文图 11/12 的验证散点拟合统计量。
%   自由斜率最小二乘：retrieved = slope*observed + intercept。
%   固定斜率为 1：retrieved = observed + bias。

    observed = double(observed(:));
    retrieved = double(retrieved(:));
    if numel(observed) ~= numel(retrieved)
        error('fit_retrieval_2017:Length', '观测和反演结果长度不一致。');
    end
    valid = isfinite(observed) & isfinite(retrieved);
    x = observed(valid);
    y = retrieved(valid);
    fit = struct('n', numel(x), 'slope', NaN, 'intercept', NaN, ...
        'correlation', NaN, 'rSquared', NaN, 'bias', NaN, 'rmse', NaN);
    if numel(x) < 2, return; end

    fit.bias = mean(y - x);  % 固定斜率为 1 时的截距。
    fit.rmse = sqrt(mean((y - x).^2));
    if std(x) > 0
        coefficients = polyfit(x, y, 1);
        fit.slope = coefficients(1);
        fit.intercept = coefficients(2);
    end
    if std(x) > 0 && std(y) > 0
        c = corrcoef(x, y);
        fit.correlation = c(1,2);
        fit.rSquared = fit.correlation^2;
    end
end
