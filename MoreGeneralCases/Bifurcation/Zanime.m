% realtime_euler_scan_style_anim_refTheta32.m
% 以 Theta32 为参考系显示动画（仅改变绘制坐标；动力学不变）

clear; clc; close all;

%% ===== 参数（与第二个程序保持一致的默认值）=====
gamma1 = 0.9;  gamma2 = 0.1;  gamma3 = 0;
Delta  = 1.0;
dt     = 1e-6;
Tmax   = 100;
K      = 50;
omega0 = 30;

plotEvery = 50;

%% ===== 读取 q(|z|) 数据（保持第二个程序的插值方式：不截断，允许样条外推）=====
load('R1_910.mat');   R1_data = solutions(:);
load('R3_910.mat');   R3_data = solutions3(:);
qfun = @(r) interp1(R3_data, R1_data, r, 'spline');   % 无 [0,1] 截断

%% ===== C(theta) =====
C = @(th) gamma1*exp(1i*(th/3)) ...
        + gamma2*exp(1i*(th/3 + 2*pi/3)) ...
        + gamma3*exp(1i*(th/3 + 4*pi/3));

%% ===== 初值（与第二个程序一致）=====
R31 = 1;   Theta31 = 0;
R32 = 1;    Theta32 = 0;

% 先用初值计算 z1, R1, Theta1
z11 = qfun(R31) * (gamma1*exp(1i*Theta31/3) + gamma2*exp(1i*(Theta31/3+2*pi/3)) + gamma3*exp(1i*(Theta31/3+4*pi/3)));
z12 = qfun(R32) * (gamma1*exp(1i*Theta32/3) + gamma2*exp(1i*(Theta32/3+2*pi/3)) + gamma3*exp(1i*(Theta32/3+4*pi/3)));
z1  = 0.5*(z11 + z12);
R1  = abs(z1);
Theta1 = angle(z1);

%% ===== 画布（第一个程序风格）=====
figure('Color','w'); hold on;
th = linspace(0,2*pi,400);
plot(cos(th), sin(th), 'k-', 'LineWidth',1.2);      % 单位圆
axis equal; axis([-1 1 -1 1]); grid on; box on;
xlabel('Re(z)'); ylabel('Im(z)');

% === 标题注明参考系：Theta32 ===
title(sprintf('Explicit Euler (2nd dynamics)  [frame: \\Theta_{32}]  K=%.2f, \\omega_0=%.2f, \\Delta=%.1f, dt=%.1e', ...
      K, omega0, Delta, dt));

traj_plus  = animatedline('Color',[0.85 0.1 0.1],'LineWidth',1.6);
traj_minus = animatedline('Color',[0.1 0.1 0.9],'LineWidth',1.6);

pt_plus  = plot(NaN,NaN,'ro','MarkerFaceColor',[1 0.6 0.6]);
pt_minus = plot(NaN,NaN,'bo','MarkerFaceColor',[0.6 0.6 1]);

% === 图例注明为 Theta32 参考系 ===
legend({'|z|=1','z_{3+} in frame of \Theta_{32}','z_{3-} in frame of \Theta_{32}'}, 'Location','southoutside');

%% ===== 显式欧拉 + 实时动画（动力学不变）=====
t = 0; k = 0;
while t < Tmax
    k = k + 1;

    % ---- 用当前(R31,Theta31,R32,Theta32)先算 z1 -> R1,Theta1 ----
    z11 = qfun(R31) * (gamma1*exp(1i*Theta31/3) + gamma2*exp(1i*(Theta31/3+2*pi/3)) + gamma3*exp(1i*(Theta31/3+4*pi/3)));
    z12 = qfun(R32) * (gamma1*exp(1i*Theta32/3) + gamma2*exp(1i*(Theta32/3+2*pi/3)) + gamma3*exp(1i*(Theta32/3+4*pi/3)));
    z1  = 0.5*(z11 + z12);
    R1  = abs(z1);
    Theta1 = angle(z1);

    % ---- ODE（无任何 R->0 保护、无 R 投影）----
    dR31     =  3*K/2 * R1^3 * cos(3*Theta1 - Theta31) * (1 - R31^2) - 3*R31*Delta;
    dTheta31 =  3*K/2 * R1^3 * sin(3*Theta1 - Theta31) * (R31 + 1/R31) - 3*omega0;
    dR32     =  3*K/2 * R1^3 * cos(3*Theta1 - Theta32) * (1 - R32^2) - 3*R32*Delta;
    dTheta32 =  3*K/2 * R1^3 * sin(3*Theta1 - Theta32) * (R32 + 1/R32) + 3*omega0;

    % ---- 显式欧拉步进 ----
    R31     = R31     + dt*dR31;
    Theta31 = Theta31 + dt*dTheta31;
    R32     = R32     + dt*dR32;
    Theta32 = Theta32 + dt*dTheta32;
    t = t + dt;

    % ---- 绘图（在 Theta32 参考系：所有角度减去 Theta32）----
    th1_rel = (Theta31 - Theta32)/3;   % z_{3+} 的相对角
    th2_rel = 0;                   % z_{3-} 在自身参考系恒为 0

    x1 = R31*cos(th1_rel);  y1 = R31*sin(th1_rel);
    x2 = R32*cos(th2_rel);  y2 = R32*sin(th2_rel);  % 即 (R32, 0)

    if mod(k, plotEvery) == 0
        addpoints(traj_plus,  x1, y1);
        addpoints(traj_minus, x2, y2);
        set(pt_plus,  'XData', x1, 'YData', y1);
        set(pt_minus, 'XData', x2, 'YData', y2);
        drawnow limitrate;
    end

    % ---- 数值异常时停止 ----
    if ~isfinite(R31) || ~isfinite(R32) || ~isfinite(Theta31) || ~isfinite(Theta32)
        warning('数值发散/异常，停止。'); break;
    end
end
