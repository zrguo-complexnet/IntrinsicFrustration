% find_saddle_equilibria.m
% 目的：在给定 K 的条件下，数值查找鞍点，并做稳定性分类。
% 依赖：R1_910.mat, R3_910.mat（与你现有代码一致）

% clear; clc; close all;

%% ===== 参数设置 =====
gamma1 = 0.9;  gamma2 = 0.1;  gamma3 = 0.0;
Delta  = 1.0;
omega0 = 70;

% 选择待检测的 K（建议在 SN 临界值附近的两侧都试一试）
K = 177;

% 求根/差分数值参数
tol_root    = 1e-10;         % 终止阈值 ||F||_2
max_newton  = 40;            % 牛顿最大迭代
fd_step_3d  = [1e-6;1e-6;1e-6]; % 数值雅可比步长
alpha_list  = [1e-3, 5e-3, 1e-2, 5e-2]; % 沿特征向量的试探步长
n_multistart= 80;            % 兜底多启动个数

% 初值（先求一个稳定不动点）
x_init = [0.9175; 0.9835; -18.684460444436220-2.286780062262495];    % [R31; R32; phi] in Theta32 frame, no wrap

%% ===== 数据与向量场 =====
load('R1_910.mat');   R1_data = solutions(:);
load('R3_910.mat');   R3_data = solutions3(:);
qfun = @(r) interp1(R3_data, R1_data, r, 'spline');   % 允许外推

params.gamma = [gamma1,gamma2,gamma3];
params.K = K; params.Delta = Delta; params.omega0 = omega0; params.qfun = qfun;

F = @(x) rhs_ref_frame(x, params);           % x=[R31; R32; phi]
Jnum = @(x) jacobian_fd_3d(F, x, fd_step_3d);

%% ===== 第一步：先找到一个稳定不动点（作为锚点）=====
[x_stb, ok_stb] = solve_equilibrium(F, Jnum, x_init, tol_root, max_newton);
if ~ok_stb
    % fsolve 兜底
    try
        opts = optimoptions('fsolve','Display','off','FunctionTolerance',tol_root,...
                            'StepTolerance',1e-12,'OptimalityTolerance',1e-12);
        x_stb = fsolve(@(y) F(y), x_init, opts);
        ok_stb = true;
    catch
        ok_stb = false;
    end
end
assert(ok_stb, '未能收敛到稳定不动点，请调整 x_init 或 K。');

% 稳定性分类
J_stb = Jnum(x_stb);
lam_stb = eig(J_stb);
fprintf('稳定不动点 x_stb = [%.6g, %.6g, %.6g]^T\n', x_stb(1), x_stb(2), x_stb(3));
fprintf('eig(J_stb) = %s\n', sprintf('[%.3e%+.3ei] ', real(lam_stb(1)), imag(lam_stb(1))));
stb_type = classify_eig(lam_stb);
fprintf('类型(稳定/鞍/源/临界)：%s\n\n', stb_type);

%% ===== 第二步：沿“最接近 0 的特征值”的特征向量出射，搜鞍点 =====
[~, idx_min] = min(abs(real(lam_stb)));   % 选取最靠近0的特征值
[v_min, ~]   = eigvec_of(J_stb, idx_min);

cands = [];  % 存放候选不动点
for sgn = [-1, 1]
    for a = alpha_list
        x0 = x_stb + sgn * a * real(v_min);
        [x_eq, ok] = solve_equilibrium(F, Jnum, x0, tol_root, max_newton);
        if ok
            cands = [cands, x_eq]; %#ok<AGROW>
        end
    end
end

% 去重（合并近似相同的解）
eqs = unique_points(cands, 1e-6);

%% ===== 如果特征向量射线法没找到新点，用多启动兜底 =====
if isempty(eqs)
    rng(0);
    cands2 = [];
    % 在一个合理盒子内撒点（可按你的系统范围调整）
    R31_box = [0.1, 1.5];
    R32_box = [0.1, 1.5];
    phi_box = [-20, 20];   % 不 wrap，允许大范围
    for k = 1:n_multistart
        x0 = [ ...
            R31_box(1) + rand*(R31_box(2)-R31_box(1));
            R32_box(1) + rand*(R32_box(2)-R32_box(1));
            phi_box(1) + rand*(phi_box(2)-phi_box(1)) ];
        [x_eq, ok] = solve_equilibrium(F, Jnum, x0, tol_root, max_newton);
        if ok
            cands2 = [cands2, x_eq]; %#ok<AGROW>
        end
    end
    eqs = unique_points(cands2, 1e-6);
end

%% ===== 分类所有找到的不动点，筛出鞍点 =====
saddles = [];
fprintf('=== 发现的不动点列表（含稳定性） ===\n');
if isempty(eqs)
    fprintf('未额外找到不动点。\n');
else
    for i = 1:size(eqs,2)
        xi = eqs(:,i);
        Ji = Jnum(xi);
        li = eig(Ji);
        tname = classify_eig(li);
        fprintf('#%d: x=[%.6g, %.6g, %.6g], type=%s, eig=%s\n', ...
            i, xi(1), xi(2), xi(3), tname, sprintf('[%.3e%+.3ei] ', real(li(1)), imag(li(1))));
        if strcmp(tname,'saddle')
            saddles = [saddles, xi]; %#ok<AGROW>
        end
    end
end

if isempty(saddles)
    fprintf('\n未找到鞍点；可尝试：靠近临界K，缩小/扩大 alpha_list；或使用更密集多启动。\n');
else
    fprintf('\n找到 %d 个鞍点。\n', size(saddles,2));
end

%% =====（可选）简单相图标注（Theta32 帧）=====
if ~isempty(saddles)
    figure('Color','w'); hold on; axis equal; grid on; box on;
    th = linspace(0,2*pi,400); plot(cos(th),sin(th),'k-'); % 单位圆
    % 画稳定点（z_{3+}*, z_{3-}*）
    plot(x_stb(1)*cos(x_stb(3)), x_stb(1)*sin(x_stb(3)), 'go','MarkerFaceColor',[0.6 1 0.6],'DisplayName','stable z_{3+}^*');
    plot(x_stb(2), 0, 'g^','MarkerFaceColor',[0.6 1 0.6],'DisplayName','stable z_{3-}^*');
    % 鞍点（红菱形）
    for i=1:size(saddles,2)
        xs = saddles(:,i);
        plot(xs(1)*cos(xs(3)), xs(1)*sin(xs(3)), 'rd','MarkerFaceColor',[1 0.6 0.6],'DisplayName','saddle z_{3+}^*');
        plot(xs(2), 0, 'rs','MarkerFaceColor',[1 0.6 0.6],'DisplayName','saddle z_{3-}^*');
    end
    xlabel('Re(z)'); ylabel('Im(z)');
    title(sprintf('Equilibria in \\Theta_{32} frame (K=%.3f)', K));
    legend('Location','best'); 
end


%% ===== 只对鞍点：边积分边绘图 (T=1, dt=1e-4) =====
if ~isempty(saddles)
    % 画底图
    figure('Color','w'); hold on; axis equal; grid on; box on;
    th = linspace(0,2*pi,400); plot(cos(th),sin(th),'k-'); % 单位圆
    xlabel('Re(z)'); ylabel('Im(z)');
    title(sprintf('Live Euler Manifolds at Saddles (T=1, dt=1e-4), K=%.3f', K));

    % 标注参考稳定点
    plot(x_stb(1)*cos(x_stb(3)), x_stb(1)*sin(x_stb(3)), 'go','MarkerFaceColor',[0.6 1 0.6],'DisplayName','stable z_{3+}^*');
    plot(x_stb(2), 0, 'g^','MarkerFaceColor',[0.6 1 0.6],'DisplayName','stable z_{3-}^*');

    % 参数
    dt = 1e-4;         % 步长
    T  = 1.0;          % 总时间
    nsteps = round(T/dt);

    eps0   = 1e-4;     % 微扰
    Rmin   = 1e-6;     % 下限
    Rmax   = 1.0;      % 上限
    phiMax = 40;       % |phi| 上限

    colU = [0.85 0.1 0.1];  % 不稳定（红）
    colS = [0.1 0.2 0.85];  % 稳定（蓝）

    % 遍历鞍点（你的 saddles 已经只含 saddle 点）
    for i = 1:size(saddles,2)
        xs = saddles(:,i);

        % 鞍点投影标注
        plot(xs(1)*cos(xs(3)), xs(1)*sin(xs(3)), 'kd','MarkerFaceColor',[1 0.8 0.8],'DisplayName','saddle z_{3+}^*');
        plot(xs(2), 0, 'ks','MarkerFaceColor',[1 0.8 0.8],'DisplayName','saddle z_{3-}^*');

        % 雅可比 & 特征方向
        J  = Jnum(xs);
        [V,D] = eig(J);
        lams = diag(D);
        idx_pos = find(real(lams) >  1e-10);
        idx_neg = find(real(lams) < -1e-10);

        v_u = []; v_s = [];
        if ~isempty(idx_pos)
            v_u = real(V(:,idx_pos(1))); v_u = v_u/max(1e-12,norm(v_u));
        end
        if ~isempty(idx_neg)
            v_s = real(V(:,idx_neg(1))); v_s = v_s/max(1e-12,norm(v_s));
        end

        % 不稳定流形：沿 v_u 的 ± 两支，正向（dt>0），实时绘图
        if ~isempty(v_u)
            for sgn = [-1, 1]
                x0 = xs + sgn*eps0*v_u;
                live_euler_and_plot(F, x0, dt, nsteps, Rmin, Rmax, phiMax, colU);
            end
        end

        % 稳定流形：沿 v_s 的 ± 两支，反向（-dt），实时绘图
        if ~isempty(v_s)
            for sgn = [-1, 1]
                x0 = xs + sgn*eps0*v_s;
                live_euler_and_plot(F, x0, -dt, nsteps, Rmin, Rmax, phiMax, colS);
            end
        end
    end
    legend('Location','best');
end

function live_euler_and_plot(F, x0, dt, nsteps, Rmin, Rmax, phiMax, col)
% 边用欧拉法积分边在投影平面画线；包含 z_{3+} 圆盘投影与 z_{3-} 实轴轨迹
    x = x0;
    % 初始化两条 animatedline：盘内曲线 和 实轴曲线
    h1 = animatedline('Color',col,'LineWidth',1.4,'MaximumNumPoints',nsteps,'DisplayName','');
    h2 = animatedline('Color',col,'LineWidth',1.0,'MaximumNumPoints',nsteps,'HandleVisibility','off');

    % 先加起点
    [xp,yp,xm] = proj_one(x);
    addpoints(h1, xp, yp);
    addpoints(h2, xm, 0);

    for k = 1:nsteps
        fx = F(x);
        x  = x + dt * fx;

        % 健壮性护栏
        if ~all(isfinite(x)), break; end
        if x(1) < Rmin || x(2) < Rmin, break; end
        if x(1) > Rmax || x(2) > Rmax, break; end
        if abs(x(3)) > phiMax,         break; end
        

        [xp,yp,xm] = proj_one(x);
        addpoints(h1, xp, yp);
        addpoints(h2, xm, 0);

        % 实时刷新（limitrate 不卡顿）
        drawnow limitrate
        % pause(0.001);
    end
end

function [xp,yp,xm] = proj_one(x)
% (R31,R32,phi) -> (Re(z_{3+}), Im(z_{3+})) 与 z_{3-} 的实轴位置
    R31 = x(1); R32 = x(2); phi = x(3);
    xp = R31*cos(phi/3);
    yp = R31*sin(phi/3);
    xm = R32;
end


%% ====== 在鞍点处用欧拉法绘制稳定/不稳定流形 ======
function plot_manifolds_at_saddle(xs, F, Jfun, opts, style)
% xs:   该鞍点坐标 [R31; R32; phi]
% F:    右端项句柄
% Jfun: 数值雅可比
% opts: 结构体，含步长/步数/半径等
% style: 结构体，控制绘图开关

    arguments
        xs (3,1) double
        F function_handle
        Jfun function_handle
        opts.h_forward (1,1) double = 1e-3      % 正向欧拉步长
        opts.h_backward (1,1) double = -1e-3     % 反向欧拉步长（负）
        opts.nsteps_u (1,1) double = 8000        % 不稳定曲线步数上限
        opts.nsteps_s (1,1) double = 4000        % 稳定曲线步数上限（每个方向）
        opts.eps0 (1,1) double = 1e-4            % 初始扰动大小
        opts.n_directions_s (1,1) double = 24    % 稳定面上的方向数
        opts.Rmin (1,1) double = 1e-6            % R 的下限，避免除零
        opts.Rmax (1,1) double = 3               % R 的上限，避免发散
        opts.phi_max (1,1) double = 40           % |phi| 上限
        opts.keep_every (1,1) double = 5         % 每隔多少步保留一个点，用于降采样
        style.plot3D (1,1) logical = false
        style.plotProjected (1,1) logical = true
        style.colorU (1,3) double = [0.85 0.1 0.1]
        style.colorS (1,3) double = [0.1 0.2 0.85]
    end

    % = 特征分解 =
    J = Jfun(xs);
    [V,D] = eig(J);
    lam   = diag(D);

    % 分类
    is_pos = real(lam) >  1e-10;
    is_neg = real(lam) < -1e-10;
    is_zero= ~(is_pos | is_neg);

    Vu = V(:,is_pos);        % 不稳定子空间
    Vs = V(:,is_neg | is_zero); % 稳定/临界子空间（临界近似并入稳定）

    % 规范化向量
    normalize = @(v) v / max(1e-12, norm(real(v)));
    if ~isempty(Vu),  Vu = real(Vu);  Vu(:,1) = normalize(Vu(:,1)); end
    if ~isempty(Vs),  Vs = real(Vs);  Vs = orth(Vs);               end  % 用 orth 得到实正交基

    % ========== 不稳定流形（通常 1 维，两条分支） ==========
    if ~isempty(Vu)
        vu = Vu(:,1); % 若>1维，可逐个分支循环
        for sgn = [-1, 1]
            x0 = xs + sgn * opts.eps0 * vu;
            Xu = euler_traj(F, x0, opts.h_forward, opts.nsteps_u, opts);
            draw_traj(Xu, style.colorU, style);
        end
    end

    % ========== 稳定流形（通常 2 维：在平面内取多方向，反向积分） ==========
    if size(Vs,2)>=1
        if size(Vs,2)==1
            % 只有 1 维稳定：两条分支（反向积分）
            vs1 = Vs(:,1);
            vs1 = normalize(vs1);
            for sgn = [-1, 1]
                x0 = xs + sgn * opts.eps0 * vs1;
                Xs = euler_traj(F, x0, opts.h_backward, opts.nsteps_s, opts);
                draw_traj(Xs, style.colorS, style);
            end
        else
            % 2 维稳定：取圆周多方向
            vs1 = Vs(:,1); vs2 = Vs(:,2);
            vs1 = normalize(vs1); vs2 = normalize(vs2);
            thetas = linspace(0, 2*pi, opts.n_directions_s+1); thetas(end)=[];
            for th = thetas
                dir = cos(th)*vs1 + sin(th)*vs2;
                x0  = xs + opts.eps0 * dir;
                Xs  = euler_traj(F, x0, opts.h_backward, opts.nsteps_s, opts);
                draw_traj(Xs, style.colorS, style);
            end
        end
    end

    % 把鞍点也标出来
    if style.plotProjected
        [xp,yp,xm] = proj_point(xs);
        plot(xp,yp,'kd','MarkerFaceColor',[1 0.8 0.8],'DisplayName','saddle z_{3+}^*');
        plot(xm,0 ,'ks','MarkerFaceColor',[1 0.8 0.8],'DisplayName','saddle z_{3-}^*');
    end
    if style.plot3D
        plot3(xs(1), xs(2), xs(3), 'kd','MarkerFaceColor',[1 0.8 0.8]);
    end
end

function X = euler_traj(F, x0, h, nsteps, opts)
% 简单欧拉积分（可正可负 h）
    X = zeros(3, floor(nsteps/opts.keep_every)+1);
    x = x0; t = 0; kkeep = 1; X(:,kkeep) = x;

    for k = 1:nsteps
        fx = F(x);

        % 欧拉一步
        x = x + h * fx;
        t = t + h;

        % 物理/数值限制，避免发散和除零
        if ~all(isfinite(x)), break; end
        if x(1) < opts.Rmin || x(2) < opts.Rmin, break; end
        if x(1) > opts.Rmax || x(2) > opts.Rmax, break; end
        if abs(x(3)) > opts.phi_max, break; end

        % 降采样记录
        if mod(k, opts.keep_every)==0
            kkeep = kkeep + 1;
            X(:,kkeep) = x;
        end
    end
    X = X(:,1:kkeep);
end

function draw_traj(X, col, style)
    if style.plotProjected
        [xp,yp,xm] = proj_points(X);
        plot(xp,yp,'-','Color',col,'LineWidth',1.2);
        plot(xm,0* xm,'-','Color',col,'LineWidth',1.2);
    end
    if style.plot3D
        plot3(X(1,:), X(2,:), X(3,:), '-', 'Color', col, 'LineWidth', 1.2);
    end
end

function [xp,yp,xm] = proj_points(X)
% 把 (R31,R32,phi) 投到你现有图的两个坐标：
% z_{3+} = R31*e^{i*phi} => (xp,yp)
% z_{3-} = R32           => (xm,0)
    R31 = X(1,:); R32 = X(2,:); phi = X(3,:);
    xp = R31 .* cos(phi);
    yp = R31 .* sin(phi);
    xm = R32;
end

function [xp,yp,xm] = proj_point(x)
    [xp,yp,xm] = proj_points(x);
end


%% ====== 函数区 ======
function F = rhs_ref_frame(x, params)
    % x = [R31; R32; phi], 在 Theta32 帧（Theta32=0, Theta31=phi）
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

    % 动力学（与你的程序一致，不做 wrap / 不做 R->0 保护）
    dR31   =  3*K/2 * R1^3 * cos(3*th1 - th31) * (1 - R31^2) - 3*R31*Delta;
    dth31  =  3*K/2 * R1^3 * sin(3*th1 - th31) * (R31 + 1/R31) - 3*omega0;
    dR32   =  3*K/2 * R1^3 * cos(3*th1 - th32) * (1 - R32^2) - 3*R32*Delta;
    dth32  =  3*K/2 * R1^3 * sin(3*th1 - th32) * (R32 + 1/R32) + 3*omega0;

    dphi   = dth31 - dth32;
    F = [dR31; dR32; dphi];
end

function [x, ok] = solve_equilibrium(F, Jfun, x0, tol_root, maxit)
    x = x0; ok = false;
    for it=1:maxit
        Fx = F(x);
        if norm(Fx) < tol_root, ok=true; break; end
        Jx = Jfun(x);
        % 牛顿修正
        dx = - Jx \ Fx;
        x  = x + dx;
        if ~all(isfinite(x)), ok=false; return; end
    end
end

function J = jacobian_fd_3d(f, x, h)
    if numel(h)==1, h = h*ones(3,1); end
    J = zeros(3,3);
    for j=1:3
        xp = x; xm = x;
        xp(j) = xp(j) + h(j);
        xm(j) = xm(j) - h(j);
        fp = f(xp); fm = f(xm);
        J(:,j) = (fp - fm)/(2*h(j));
    end
end

function v = normalize_vec(v)
    n = norm(v); if n>0, v = v/n; end
end

function [v, D] = eigvec_of(J, idx)
    [V,D] = eig(J);
    v = V(:,idx);
    v = normalize_vec(v);
end

function T = classify_eig(lams)
    rp = real(lams);
    has_pos = any(rp > 0);
    has_neg = any(rp < 0);
    near0   = any(abs(rp) < 1e-8);
    if has_pos && has_neg
        T = 'saddle';
    elseif ~has_pos && any(rp < 0)
        T = 'stable';
    elseif has_pos && ~has_neg
        T = 'source';
    elseif near0
        T = 'critical';
    else
        T = 'unknown';
    end
end

function Xuniq = unique_points(X, tol)
    if isempty(X), Xuniq = []; return; end
    keep = true(1,size(X,2));
    for i=1:size(X,2)
        if ~keep(i), continue; end
        for j=i+1:size(X,2)
            if norm(X(:,i)-X(:,j)) < tol
                keep(j) = false;
            end
        end
    end
    Xuniq = X(:,keep);
end
