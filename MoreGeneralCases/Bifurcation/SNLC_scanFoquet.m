% 并行扫描 K：从 110 递减到 95，计算每个 K 的 Floquet 乘子（失败自动跳过）
clc; clear;

%% ===== 固定参数（与你当前脚本一致）=====
gamma1 = 0.9; %#ok<NASGU>
gamma2 = 0.1; %#ok<NASGU>
gamma3 = 0.0; %#ok<NASGU>

omega0 = 100;
Delta  = 1;
dt     = 1e-5;
Tmax   = 5;

% 初值（每个 K 相同）
R31_0 = 1.0;  Theta31_0 = 0.0;
R32_0 = 1.0;  Theta32_0 = 0.0;

% R1(R3) 数据与插值器（用 griddedInterpolant 更适合并行）
load('R1_910.mat');  % solutions -> R1_data
load('R3_910.mat');  % solutions3 -> R3_data
R1_data = solutions;
R3_data = solutions3;
R1_from_R3 = @(R3) interp1(R3_data, R1_data, R3, 'spline', 'extrap');

% Poincaré & Floquet 参数（保持不变）
section_dir = -1;
M_settle     = 6;
tol_section  = 1e-3;
perturb_var  = 'R31';
eps_perturb  = 1e-1;
verbose      = false; % 并行时关闭打印

% 扫描序列（从大到小）
K_vec = 110:-0.2:95;
nK = numel(K_vec);

% 结果容器（并行-friendly）
mu_all = nan(nK,1);
K_ok   = nan(nK,1);

%% ===== 并行池 =====
if isempty(gcp('nocreate'))
    parpool; % 默认开最大可用核
end

%% ===== 并行循环 =====
parfor idx = 1:nK
    K = K_vec(idx);
    try
        [ok, mu] = run_one_K(K, omega0, Delta, dt, Tmax, ...
                             R31_0, Theta31_0, R32_0, Theta32_0, ...
                             R1_from_R3, section_dir, M_settle, tol_section, ...
                             perturb_var, eps_perturb, verbose);
        if ok && isfinite(mu) && ~isnan(mu)
            mu_all(idx) = mu;
            K_ok(idx)   = K;
        end
    catch
        % 该 K 失败就跳过
    end
end

%% ===== 收集有效结果并绘图 =====
mask = isfinite(mu_all);
K_keep  = K_ok(mask);
mu_keep = mu_all(mask);

figure('Color','w');
plot(K_keep, mu_keep, 'o-','LineWidth',1.5,'MarkerSize',5);
set(gca,'XDir','reverse'); % 从 110 降到 95
xlabel('K'); ylabel('\mu (Floquet multiplier)');
title('Floquet multiplier vs. K (parallel)');
grid on;


%% ===================== 辅助函数 =====================

function [ok, mu] = run_one_K(K, omega0, Delta, dt, Tmax, ...
                              R31, Theta31, R32, Theta32, ...
                              R1_from_R3, section_dir, M_settle, tol_section, ...
                              perturb_var, eps_perturb, verbose)
    % 单个 K 的完整流程：判稳极限环 → 计算乘子
    ok = false; mu = NaN;

    Nmax = round(Tmax/dt);
    phi_prev = wrapToPi(Theta31 - Theta32);
    t = 0;

    R_section = [];
    state_at_section = [];
    t_at_section = [];

    % —— 主积分，寻找极限环并在截面上判稳（与原逻辑一致）——
    for k = 1:Nmax
        t = t + dt;

        [k1_R31,k1_T31,k1_R32,k1_T32] = rhs(R31,Theta31,R32,Theta32,K,omega0,Delta,R1_from_R3);
        [k2_R31,k2_T31,k2_R32,k2_T32] = rhs(R31+0.5*dt*k1_R31, Theta31+0.5*dt*k1_T31, ...
                                            R32+0.5*dt*k1_R32, Theta32+0.5*dt*k1_T32, ...
                                            K,omega0,Delta,R1_from_R3);
        [k3_R31,k3_T31,k3_R32,k3_T32] = rhs(R31+0.5*dt*k2_R31, Theta31+0.5*dt*k2_T31, ...
                                            R32+0.5*dt*k2_R32, Theta32+0.5*dt*k2_T32, ...
                                            K,omega0,Delta,R1_from_R3);
        [k4_R31,k4_T31,k4_R32,k4_T32] = rhs(R31+dt*k3_R31, Theta31+dt*k3_T31, ...
                                            R32+dt*k3_R32, Theta32+dt*k3_T32, ...
                                            K,omega0,Delta,R1_from_R3);

        R31     = R31     + dt/6*(k1_R31 + 2*k2_R31 + 2*k3_R31 + k4_R31);
        Theta31 = Theta31 + dt/6*(k1_T31 + 2*k2_T31 + 2*k3_T31 + k4_T31);
        R32     = R32     + dt/6*(k1_R32 + 2*k2_R32 + 2*k3_R32 + k4_R32);
        Theta32 = Theta32 + dt/6*(k1_T32 + 2*k2_T32 + 2*k3_T32 + k4_T32);

        phi = wrapToPi(Theta31 - Theta32);
        crossed = (section_dir>0 && (phi_prev < 0 && phi >= 0 && cos(phi/3)>0)) || ...
                  (section_dir<0 && (phi_prev > 0 && phi <= 0 && cos(phi/3)>0));
        if crossed
            alpha = abs(phi_prev) / (abs(phi_prev) + abs(phi) + eps);
            R31_sec     = R31     - alpha * (R31 - (R31 - dt/6*(k1_R31 + 2*k2_R31 + 2*k3_R31 + k4_R31)));
            Theta31_sec = Theta31 - alpha * (Theta31 - (Theta31 - dt/6*(k1_T31 + 2*k2_T31 + 2*k3_T31 + k4_T31)));
            R32_sec     = R32     - alpha * (R32 - (R32 - dt/6*(k1_R32 + 2*k2_R32 + 2*k3_R32 + k4_R32)));
            Theta32_sec = Theta32 - alpha * (Theta32 - (Theta32 - dt/6*(k1_T32 + 2*k2_T32 + 2*k3_T32 + k4_T32)));
            t_sec       = t - alpha*dt;

            switch perturb_var
                case 'R31', R_on_sec = R31_sec;
                case 'R32', R_on_sec = R32_sec;
                otherwise, error('perturb_var must be ''R31'' or ''R32''.');
            end
            R_section(end+1,1) = R_on_sec; %#ok<AGROW>
            state_at_section(end+1,:) = [R31_sec,Theta31_sec,R32_sec,Theta32_sec]; %#ok<AGROW>
            t_at_section(end+1,1) = t_sec; %#ok<AGROW>

            if numel(R_section) >= M_settle+1
                recent = R_section(end-M_settle:end);
                if max(abs(diff(recent))) < tol_section
                    % 判稳
                    R_star = R_section(end);
                    x_star = state_at_section(end,:);

                    % 计算 Floquet 乘子（中心差分，保持原逻辑）
                    mu = computeFloquet(x_star, perturb_var, eps_perturb, ...
                                        K,omega0,Delta,dt,section_dir,R1_from_R3, ...
                                        [],[],[],[],Inf,verbose);
                    if isfinite(mu) && ~isnan(mu)
                        ok = true;
                    end
                    return;
                end
            end
        end
        phi_prev = phi;
    end
end

% ===== 下面三个函数与原脚本一致（去除动画参数并保持兼容）=====
function [dR31,dT31,dR32,dT32] = rhs(R31,Theta31,R32,Theta32,K,omega0,Delta,R1_from_R3)
    gamma1 = 0.9; gamma2 = 0.1; gamma3 = 0.0;
    z11 = R1_from_R3(R31) * ( ...
        gamma1*exp(1i*(Theta31/3)) + ...
        gamma2*exp(1i*(Theta31/3 + 2*pi/3)) + ...
        gamma3*exp(1i*(Theta31/3 + 4*pi/3)) );
    z12 = R1_from_R3(R32) * ( ...
        gamma1*exp(1i*(Theta32/3)) + ...
        gamma2*exp(1i*(Theta32/3 + 2*pi/3)) + ...
        gamma3*exp(1i*(Theta32/3 + 4*pi/3)) );
    z1 = 0.5*(z11 + z12);
    R1 = abs(z1); %#ok<NASGU>
    Theta1 = angle(z1);

    R31 = max(R31, 1e-12);
    R32 = max(R32, 1e-12);

    dR31 = (3*K/2) * abs(z1)^3 * cos(3*Theta1 - Theta31) * (1 - R31^2) - 3*Delta*R31;
    dT31 = (3*K/2) * abs(z1)^3 * sin(3*Theta1 - Theta31) * (R31 + 1/R31) - 3*omega0;

    dR32 = (3*K/2) * abs(z1)^3 * cos(3*Theta1 - Theta32) * (1 - R32^2) - 3*Delta*R32;
    dT32 = (3*K/2) * abs(z1)^3 * sin(3*Theta1 - Theta32) * (R32 + 1/R32) + 3*omega0;
end

function mu = computeFloquet(x_star, perturb_var, eps0, ...
                             K,omega0,Delta,dt,section_dir,R1_from_R3, ...
                             ~,~,~,~,~,verbose)
    eps = eps0;
    x_plus  = x_star; x_minus = x_star;
    switch perturb_var
        case 'R31'
            x_plus(1)  = max(1e-12, x_star(1)  + eps);
            x_minus(1) = max(1e-12, x_star(1)  - eps);
        case 'R32'
            x_plus(3)  = max(1e-12, x_star(3)  + eps);
            x_minus(3) = max(1e-12, x_star(3)  - eps);
        otherwise
            error('perturb_var must be ''R31'' or ''R32''.');
    end

    Rp_next = advance_to_next_section(x_plus(1),x_plus(2),x_plus(3),x_plus(4), ...
        K,omega0,Delta,dt,section_dir,R1_from_R3);
    Rm_next = advance_to_next_section(x_minus(1),x_minus(2),x_minus(3),x_minus(4), ...
        K,omega0,Delta,dt,section_dir,R1_from_R3);

    switch perturb_var
        case 'R31', mu = (Rp_next(1) - Rm_next(1)) / (2*eps);
        case 'R32', mu = (Rp_next(3) - Rm_next(3)) / (2*eps);
    end

    if verbose, fprintf('  (central diff) collected.\n'); end
end

function x_next = advance_to_next_section(R31,Theta31,R32,Theta32, ...
        K,omega0,Delta,dt,section_dir,R1_from_R3)
    % 离开切面一步
    [k1_R31,k1_T31,k1_R32,k1_T32] = rhs(R31,Theta31,R32,Theta32,K,omega0,Delta,R1_from_R3);
    [k2_R31,k2_T31,k2_R32,k2_T32] = rhs(R31+0.5*dt*k1_R31, Theta31+0.5*dt*k1_T31, ...
                                        R32+0.5*dt*k1_R32, Theta32+0.5*dt*k1_T32, ...
                                        K,omega0,Delta,R1_from_R3);
    [k3_R31,k3_T31,k3_R32,k3_T32] = rhs(R31+0.5*dt*k2_R31, Theta31+0.5*dt*k2_T31, ...
                                        R32+0.5*dt*k2_R32, Theta32+0.5*dt*k2_T32, ...
                                        K,omega0,Delta,R1_from_R3);
    [k4_R31,k4_T31,k4_R32,k4_T32] = rhs(R31+dt*k3_R31, Theta31+dt*k3_T31, ...
                                        R32+dt*k3_R32, Theta32+dt*k3_T32, ...
                                        K,omega0,Delta,R1_from_R3);
    R31     = R31     + dt/6*(k1_R31 + 2*k2_R31 + 2*k3_R31 + k4_R31);
    Theta31 = Theta31 + dt/6*(k1_T31 + 2*k2_T31 + 2*k3_T31 + k4_T31);
    R32     = R32     + dt/6*(k1_R32 + 2*k2_R32 + 2*k3_R32 + k4_R32);
    Theta32 = Theta32 + dt/6*(k1_T32 + 2*k2_T32 + 2*k3_T32 + k4_T32);

    phi_prev = wrapToPi(Theta31 - Theta32);
    while true
        [k1_R31,k1_T31,k1_R32,k1_T32] = rhs(R31,Theta31,R32,Theta32,K,omega0,Delta,R1_from_R3);
        [k2_R31,k2_T31,k2_R32,k2_T32] = rhs(R31+0.5*dt*k1_R31, Theta31+0.5*dt*k1_T31, ...
                                            R32+0.5*dt*k1_R32, Theta32+0.5*dt*k1_T32, ...
                                            K,omega0,Delta,R1_from_R3);
        [k3_R31,k3_T31,k3_R32,k3_T32] = rhs(R31+0.5*dt*k2_R31, Theta31+0.5*dt*k2_T31, ...
                                            R32+0.5*dt*k2_R32, Theta32+0.5*dt*k2_T32, ...
                                            K,omega0,Delta,R1_from_R3);
        [k4_R31,k4_T31,k4_R32,k4_T32] = rhs(R31+dt*k3_R31, Theta31+dt*k3_T31, ...
                                            R32+dt*k3_R32, Theta32+dt*k3_T32, ...
                                            K,omega0,Delta,R1_from_R3);

        R31_new     = R31     + dt/6*(k1_R31 + 2*k2_R31 + 2*k3_R31 + k4_R31);
        Theta31_new = Theta31 + dt/6*(k1_T31 + 2*k2_T31 + 2*k3_T31 + k4_T31);
        R32_new     = R32     + dt/6*(k1_R32 + 2*k2_R32 + 2*k3_R32 + k4_R32);
        Theta32_new = Theta32 + dt/6*(k1_T32 + 2*k2_T32 + 2*k3_T32 + k4_T32);

        phi = wrapToPi(Theta31_new - Theta32_new);
        crossed = (section_dir>0 && (phi_prev < 0 && phi >= 0 && cos(phi/3)>0)) || ...
                  (section_dir<0 && (phi_prev > 0 && phi <= 0 && cos(phi/3)>0));
        if crossed
            alpha = abs(phi_prev) / (abs(phi_prev) + abs(phi) + eps);
            R31_sec     = R31_new     - alpha * (R31_new     - R31);
            Theta31_sec = Theta31_new - alpha * (Theta31_new - Theta31);
            R32_sec     = R32_new     - alpha * (R32_new     - R32);
            Theta32_sec = Theta32_new - alpha * (Theta32_new - Theta32);
            x_next = [R31_sec,Theta31_sec,R32_sec,Theta32_sec];
            return;
        end

        R31=R31_new; Theta31=Theta31_new; R32=R32_new; Theta32=Theta32_new; phi_prev=phi;
    end
end

function a = wrapToPi(a)
    % 保持你当前的定义（区间宽度为 2π 的对称包装）
    a = mod(a + 3*pi, 6*pi) - 3*pi;
end
