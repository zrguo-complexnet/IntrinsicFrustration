% fit_floquet_sqrt.m
% 对 SNLC 附近 Floquet 乘子使用平方根模型 lambda = 1 - A*sqrt(mu - mu_c) 进行拟合并作图
% MATLAB R2018a+ 兼容。优先使用 lsqcurvefit（Optimization Toolbox），若无则回退到 fminsearch。

%% ========== 数据读取 ==========
% 方式A：从CSV读取（两列：mu, lambda）。可带表头。
% dataFile = 'floquet_data.csv';
% mu = []; lambdaVals = [];
% 
% if exist(dataFile,'file')
%     T = readtable(dataFile);
%     if width(T) < 2
%         error('CSV 至少需要两列：第一列 mu，第二列 lambda。');
%     end
%     mu = T{:,1};
%     lambdaVals = T{:,2};
% end

% 方式B：【可选：直接粘贴数据】如果你不想读CSV，把下面两行解除注释并填入向量：
mu = K_keep;
lambdaVals = mu_keep;

if isempty(mu)
    error(['未找到数据。请提供 CSV 文件 "' dataFile '"，或在脚本中直接赋值 mu 与 lambdaVals。']);
end

mu = mu(:); lambdaVals = lambdaVals(:);
[mu, idx] = sort(mu);
lambdaVals = lambdaVals(idx);

%% ===== 2) 固定 Kc =====
Kc = 105.9;   % <<<=== 这里改成你已知的 Kc 数值 ===>>>

fprintf('使用固定 Kc = %.6g 进行拟合。\n', Kc);

%% ===== 3) 拟合 A =====
% 模型: λ = 1 - A * sqrt(μ - Kc)
model = @(A, m) 1 - A * sqrt(max(m - Kc, 0));

% 只拟合 A
opts = optimoptions('lsqcurvefit','Display','off');
A0 = max((1 - min(lambdaVals)) / sqrt(max(max(mu) - Kc, 1e-6)), 1e-4);
lb = [0]; ub = [Inf];
[A_fit, ~] = lsqcurvefit(@(A,m) model(A,m), A0, mu, lambdaVals, lb, ub, opts);

%% ===== 4) 绘图（带延拓）=====
figure('Color','w','Position',[100 100 800 600]);
hold on; box on;

% 原始数据
scatter(mu, lambdaVals, 40, 'filled', 'DisplayName', '数据');

% 延拓比例
extend_ratio = 0.1;
muRange = max(mu) - min(mu);
muFine = linspace(min(mu) - extend_ratio*muRange, ...
                  max(mu) + extend_ratio*muRange, 400);

% 拟合曲线
plot(muFine, model(A_fit, muFine), 'r-', 'LineWidth', 2, 'DisplayName', '拟合曲线');
yline(1, '--k', 'DisplayName', '\lambda = 1');

xlabel('\mu');
ylabel('\lambda(\mu)');
title(sprintf('固定 K_c = %.3g 的平方根拟合', Kc), 'FontSize', 12);
legend('Location','best');

% 参数标注
txt = sprintf('A = %.4g\nK_c = %.4g (固定)', A_fit, Kc);
xl = xlim; yl = ylim;
text(xl(1) + 0.02*(xl(2)-xl(1)), yl(1) + 0.05*(yl(2)-yl(1)), txt, ...
    'FontSize', 11, 'BackgroundColor',[0.95 0.95 0.98], 'EdgeColor',[0.8 0.8 0.9]);

saveas(gcf, 'fit_floquet_fixedKc.png');
fprintf('绘图完成，结果已保存为 fit_floquet_fixedKc.png\n');