% ================== Only plot order-parameter angular velocity (Omega) ==================
% Model: harmonic Kuramoto with mean-field m-th order coupling:
%   dtheta_i/dt = omega_i + K * R^m * sin(m*(Theta - theta_i))
%   Z = R*exp(i*Theta) = (1/N) * sum_j exp(i*theta_j)
%
% Single figure with 4 lines for m=1,2,3,4. Points with tail-avg |Z| < Rmin are NaN.
% ------------------------------------------------------------------------
clear; clc; close all;
rng(1);                        % 可复现

% ------------------------ 全局参数 ------------------------
N      = 1000;                 % 振子数
sigma  = 1.0;                  % Lorentzian 宽度
T      = 20.0;                 % 总仿真时长
dt     = 1e-3;                 % 时间步长
nSteps = round(T/dt);
tvec   = (0:nSteps-1)*dt;

Rmin   = 0.1;                  % 仅当 |Z| >= Rmin 才计入Omega

% 不同阶的 K 扫描范围（按需修改）
Kgrid = struct();
Kgrid.m1 = linspace(0, 30, 30);     % m=1 的 K 取值
Kgrid.m2 = linspace(0, 30, 30);     % m=2 的 K 取值
Kgrid.m3 = linspace(0, 30, 30);     % m=3 的 K 取值
Kgrid.m4 = linspace(0, 30, 30);     % m=4 的 K 取值（可按需调大范围）

mList  = [1 2 3 4];                % 要绘制的阶数

% 初相位：0, 2pi/3, 4pi/3 按 9:1:0
ratios = [0.9 0.1 0.0];
pts    = [0, 2*pi/3, 4*pi/3];

% 末段平均比例（用于稳态统计）
tailFrac = 0.3;                     % 最后 30% 时间
tailIdx  = round((1-tailFrac)*nSteps):nSteps;

% ------------------------ 自然频率（对称单峰 Lorentzian） ------------------------
% Lorentzian: omega = sigma * tan(pi*(U - 0.5)), mean = 0
U = rand(N,1);
omega = sigma * tan(pi*(U - 0.5));

% ------------------------ 初相位 ------------------------
theta0 = zeros(N,1);
counts = round(N*ratios);
counts(end) = N - sum(counts(1:end-1));  % 保证和为 N
idx = 1;
for g = 1:3
    if counts(g) > 0
        theta0(idx:idx+counts(g)-1) = pts(g);
        idx = idx + counts(g);
    end
end

% ================== 主循环：m=1,2,3,4（同一张图绘制） ==================
figure('Color','w'); hold on; box on; grid on;
colors  = lines(numel(mList));
markers = {'o','s','^','d'};
h = gobjects(1,numel(mList));   % 曲线句柄（用于图例）
legendNames = cell(1,numel(mList));

for im = 1:numel(mList)
    m = mList(im);
    switch m
        case 1, Kvals = Kgrid.m1;
        case 2, Kvals = Kgrid.m2;
        case 3, Kvals = Kgrid.m3;
        case 4, Kvals = Kgrid.m4;
    end

    Omega_vs_K = nan(size(Kvals));   % 预分配
    Rtail_mean = nan(size(Kvals));   % 记录尾段 |Z| 平均（可选）

    fprintf('Simulating m = %d ...\n', m);
    for ik = 1:numel(Kvals)
        K = Kvals(ik);

        % 相同步长 RK4 积分
        theta = theta0;
        Ztail = zeros(1, numel(tailIdx));   % 仅存尾段以节省内存
        tcnt  = 0;

        for it = 1:nSteps
            % 当前 Z
            Z  = mean(exp(1i*theta));
            R  = abs(Z); Th = angle(Z);

            % 右端项
            f = @(th) omega + K * (R^m) .* sin( m*(Th - th) );

            % RK4
            k1 = f(theta);

            Z2 = mean(exp(1i*(theta + 0.5*dt*k1))); R2 = abs(Z2); Th2 = angle(Z2);
            f2 = @(th) omega + K * (R2^m) .* sin( m*(Th2 - th) );
            k2 = f2(theta + 0.5*dt*k1);

            Z3 = mean(exp(1i*(theta + 0.5*dt*k2))); R3 = abs(Z3); Th3 = angle(Z3);
            f3 = @(th) omega + K * (R3^m) .* sin( m*(Th3 - th) );
            k3 = f3(theta + 0.5*dt*k2);

            Z4 = mean(exp(1i*(theta + dt*k3))); R4 = abs(Z4); Th4 = angle(Z4);
            f4 = @(th) omega + K * (R4^m) .* sin( m*(Th4 - th) );
            k4 = f4(theta + dt*k3);

            theta = theta + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
            theta = angle(exp(1i*theta));   % wrap to (-pi, pi]

            % 尾段记录 Z
            if it >= tailIdx(1)
                tcnt = tcnt + 1;
                Ztail(tcnt) = mean(exp(1i*theta));
            end
        end

        % 末段统计
        Rt = abs(Ztail);
        Rtail_mean(ik) = mean(Rt);

        if Rtail_mean(ik) >= Rmin
            ph = angle(Ztail);
            ph = unwrap(ph);
            tt = tvec(tailIdx);
            p = polyfit(tt, ph, 1);
            Omega_vs_K(ik) = p(1);     % 斜率即集体角速度
        else
            Omega_vs_K(ik) = NaN;      % |Z| 太小，视为无效/不绘制
        end
    end

    % -------- 在同一坐标轴上绘制本阶的 Omega vs K --------
    h(im) = plot(Kvals, Omega_vs_K, ...
        'LineWidth', 2, 'MarkerSize', 5, ...
        'Color', colors(im,:), 'Marker', markers{im});
    legendNames{im} = sprintf('m = %d', m);
end

% 轴标签/标题/图例
yline(0,'k--','LineWidth',1);
xlabel('K'); ylabel('\Omega  (d/dt arg(Z))');
title(sprintf('\\Omega vs K for m = 1,2,3,4  (only \\langle|Z|\\rangle_{tail} \\ge %.2f)', Rmin));
legend(h, legendNames, 'Location','best');
hold off;

% ------------------------ 备注 ------------------------
% 1) 可单独调整 Kgrid.m4 的范围；高阶常需要更大的 K 才显现频率牵引。
% 2) 如需看尾段 |Z|，可另开一幅图绘制 Rtail_mean 与 K 的关系。
