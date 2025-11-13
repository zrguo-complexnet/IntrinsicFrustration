% sweep_K_Jplus_eigs.m
% 扫描 K 从 85 到 80，按照“Theta32 参考系 + 不 wrap 角度 + 中心差分雅可比”的方法：
% 1) 在每个 K 上求 Theta32 帧的不动点 x* = [R31*, R32*, phi*]
% 2) 固定 R32 = R32*，计算 z_{3+} 子系统 (u=[R31; phi]) 的 2x2 Jacobian: J_plus
% 3) 记录 eig(J_plus) 随 K 的变化并绘图（实部）

clear; clc; 

%% ===== 固定参数（除 K 外）=====
gamma1 = 0.9;  gamma2 = 0.1;  gamma3 = 0;
Delta  = 1.0;
omega0 = 40;

% 插值数据
load('R1_910.mat');   R1_data = solutions(:);
load('R3_910.mat');   R3_data = solutions3(:);
qfun = @(r) interp1(R3_data, R1_data, r, 'spline');   % 允许外推，不截断

% 连续延拓的 K 取值
K_vals = linspace(120, 105.9, 101);   % 101点，步长约0.05
nK = numel(K_vals);

% 数值设置
tol_root = 1e-10;        % 求根容差
fd_h_R   = 1e-6;         % 差分步长（R）
fd_h_phi = 1e-6;         % 差分步长（phi）
max_iter_nr = 30;        % 退路牛顿迭代次数（无 fsolve 时）

%% ===== 共享函数句柄（参考系向量场）=====
params_template.gamma = [gamma1, gamma2, gamma3];
params_template.Delta = Delta;
params_template.omega0 = omega0;
params_template.qfun = qfun;
% K 会在循环中填入

% 初始猜测：从 K=85 开始
R31_init = 1.0;  R32_init = 1.0;  Theta31_init = 0.0;  Theta32_init = 0.0;
phi_init = Theta31_init - Theta32_init;
x_guess  = [R31_init; R32_init; phi_init];   % Theta32 帧的状态变量

% 结果存储
eig1 = nan(1, nK);   % J_plus 的特征值 λ1
eig2 = nan(1, nK);   % J_plus 的特征值 λ2
maxReal = nan(1, nK);
R31_star = nan(1, nK); R32_star = nan(1, nK); phi_star = nan(1, nK);

%% ===== 主循环：对每个 K 求不动点与 J_plus =====
for idx = 1:nK
    K = K_vals(idx);
    params = params_template; params.K = K;

    % 本 K 下的向量场（Theta32 参考系）
    f_ref = @(x) rhs_ref_frame(x, params);     % x=[R31; R32; phi], 返回 [dR31; dR32; dphi]

    % ---- 先用 fsolve；失败则用简易牛顿（无需 wrap 角度）----
    x_star = x_guess;
    success = false;
    try
        opts = optimoptions('fsolve','Display','off','FunctionTolerance',tol_root,...
                            'StepTolerance',1e-12,'OptimalityTolerance',1e-12);
        x_star = fsolve(@(y) f_ref(y), x_guess, opts);
        success = true;
    catch
        % 退路：数值雅可比 + 牛顿
        x = x_guess;
        for it = 1:max_iter_nr
            F = f_ref(x);
            if norm(F) < tol_root, break; end
            Jnum = jacobian_fd_3d(f_ref, x, [1e-6; 1e-6; 1e-6]);
            % 解线性系统 J * dx = -F
            dx = - Jnum \ F;
            x  = x + dx;
            if ~all(isfinite(x)), break; end
        end
        if norm(f_ref(x)) < 1e-6
            x_star = x;
            success = true;
        end
    end

    if ~success || ~all(isfinite(x_star))
        warning('K=%.3f 未能可靠求得不动点，跳过。', K);
        % 不更新 x_guess，留给下一个 K 使用上一次成功点
        continue
    end

    % 记录不动点
    R31s = x_star(1); R32s = x_star(2); phis = x_star(3);
    R31_star(idx) = R31s; R32_star(idx) = R32s; phi_star(idx) = phis;

    % ---- 计算 z_{3+} 子系统的 2x2 Jacobian（固定 R32=R32*）----
    g_plus = @(u) sub_rhs_plus(u, R32s, params);  % u=[R31; phi] -> [dR31; dphi]
    Jplus  = jacobian_fd_2d(g_plus, [R31s; phis], [fd_h_R; fd_h_phi]);
    lams   = eig(Jplus);

    % 为了画“实部随 K”，保存
    eig1(idx) = lams(1);
    eig2(idx) = lams(2);
    maxReal(idx) = max(real(lams));

    % 连续延拓：用当前不动点作为下一个 K 的初值
    x_guess = x_star;
end

%% ===== 绘图：特征值实部随 K 的变化 =====
figure('Color','w'); hold on; grid on; box on;
plot(K_vals, real(eig1), '-', 'LineWidth', 1.8);
plot(K_vals, real(eig2), '-', 'LineWidth', 1.8);
plot(K_vals, maxReal,  '--', 'LineWidth', 1.5);   % 方便观察最危险指数
yline(0,'k:','LineWidth',1.0);

xlabel('K'); ylabel('Real part of eigenvalues');
title('Re(\lambda) of J_{plus} vs K  (Theta_{32} frame, R_{32} fixed at equilibrium)');
legend({'Re(\lambda_1)','Re(\lambda_2)','max Re(\lambda)'}, 'Location','best');

% 如果需要，标注可能的“临界K”（max Re(lambda) 最接近 0 的点）
[~, iCrit] = min(abs(maxReal));
if isfinite(K_vals(iCrit))
    xline(K_vals(iCrit),'r:','LineWidth',1.2);
    text(K_vals(iCrit), 0, sprintf('  K_c \\approx %.3f', K_vals(iCrit)), ...
        'Color',[0.8 0 0], 'VerticalAlignment','bottom');
end

%% =====（可选）把不动点也画出来便于检查 =====
figure('Color','w');
tiledlayout(3,1,'Padding','compact','TileSpacing','compact');
nexttile; plot(K_vals, R31_star,'LineWidth',1.6); grid on; ylabel('R_{31}^*'); title('Equilibria vs K');
nexttile; plot(K_vals, R32_star,'LineWidth',1.6); grid on; ylabel('R_{32}^*');
nexttile; plot(K_vals, phi_star,'LineWidth',1.6); grid on; ylabel('\phi^*'); xlabel('K');

%% ====== 向量场与数值雅可比（与前文一致，无 wrap） ======
function F = rhs_ref_frame(x, params)
    % x = [R31; R32; phi],  在 Theta32 帧（令 Theta32=0, Theta31=phi）
    R31 = x(1); R32 = x(2); phi = x(3);
    th31 = phi; th32 = 0;

    qfun = params.qfun;
    g = params.gamma;
    C = @(th) g(1)*exp(1i*(th/3)) + g(2)*exp(1i*(th/3 + 2*pi/3)) + g(3)*exp(1i*(th/3 + 4*pi/3));

    z11 = qfun(R31) * C(th31);
    z12 = qfun(R32) * C(th32);
    z1  = 0.5*(z11 + z12);
    R1  = abs(z1); th1 = angle(z1);

    K = params.K; Delta = params.Delta; omega0 = params.omega0;

    % ODE
    dR31   =  3*K/2 * R1^3 * cos(3*th1 - th31) * (1 - R31^2) - 3*R31*Delta;
    dth31  =  3*K/2 * R1^3 * sin(3*th1 - th31) * (R31 + 1/R31) - 3*omega0;
    dR32   =  3*K/2 * R1^3 * cos(3*th1 - th32) * (1 - R32^2) - 3*R32*Delta;
    dth32  =  3*K/2 * R1^3 * sin(3*th1 - th32) * (R32 + 1/R32) + 3*omega0;

    dphi   = dth31 - dth32;    % 相对相位
    F = [dR31; dR32; dphi];
end

function G = sub_rhs_plus(u, R32fix, params)
    % z_{3+} 子系统：u=[R31; phi]，R32 固定
    R31 = u(1); phi = u(2);
    x = [R31; R32fix; phi];
    F = rhs_ref_frame(x, params);
    G = [F(1); F(3)];  % 只取 dR31 与 dphi
end

function J = jacobian_fd_2d(f, u, h)
    if numel(h)==1, h = h*ones(2,1); end
    J = zeros(2,2);
    for j=1:2
        up = u; um = u;
        up(j)=up(j)+h(j); um(j)=um(j)-h(j);
        fp = f(up); fm = f(um);
        J(:,j) = (fp - fm)/(2*h(j));
    end
end

function J = jacobian_fd_3d(f, x, h)
    if numel(h)==1, h = h*ones(3,1); end
    J = zeros(3,3);
    for j=1:3
        xp = x; xm = x;
        xp(j)=xp(j)+h(j); xm(j)=xm(j)-h(j);
        fp = f(xp); fm = f(xm);
        J(:,j) = (fp - fm)/(2*h(j));
    end
end
