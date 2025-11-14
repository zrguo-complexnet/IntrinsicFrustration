% Scan K from 120 down to 105.9 using the “Theta32 reference frame + no angle wrapping
% + centered finite-difference Jacobian” approach:
%
%   1) For each value of K, solve for the equilibrium x* = [R31*, R32*, phi*]
%      in the Theta32 rotating frame.
%   2) Fix R32 = R32* and compute the 2×2 Jacobian J_plus of the z_{3+}
%      subsystem with respect to u = [R31; phi].
%   3) Record the eigenvalues of J_plus as functions of K and plot the
%      evolution of their real parts.
%

clear; clc;

%% ===== Fixed parameters (except K) =====
gamma1 = 0.9;  gamma2 = 0.1;  gamma3 = 0;
Delta  = 1.0;
omega0 = 40;

% Load interpolation data
load('R1_910.mat');   R1_data = solutions(:);
load('R3_910.mat');   R3_data = solutions3(:);
qfun = @(r) interp1(R3_data, R1_data, r, 'spline');   % Allow extrapolation, no truncation

% K values for continuation
K_vals = linspace(120, 105.9, 101);   % 101 sample points, ≈0.05 step
nK = numel(K_vals);

% Numerical settings
tol_root = 1e-10;        % Root-finding tolerance
fd_h_R   = 1e-6;         % Finite difference step (R-direction)
fd_h_phi = 1e-6;         % Finite difference step (phi-direction)
max_iter_nr = 30;        % Newton fallback iterations (if fsolve unavailable)

%% ===== Shared parameter template (K assigned inside the loop) =====
params_template.gamma = [gamma1, gamma2, gamma3];
params_template.Delta = Delta;
params_template.omega0 = omega0;
params_template.qfun = qfun;
% params_template.K will be set in the loop

% Initial guess for K = 120
R31_init = 1.0;  R32_init = 1.0;  
Theta31_init = 0.0;  
Theta32_init = 0.0;
phi_init = Theta31_init - Theta32_init;
x_guess  = [R31_init; R32_init; phi_init];   % State in Theta32 frame

% Storage variables
eig1 = nan(1, nK);          % First eigenvalue of J_plus
eig2 = nan(1, nK);          % Second eigenvalue of J_plus
maxReal = nan(1, nK);       % Maximum real part among eigenvalues

R31_star = nan(1, nK);
R32_star = nan(1, nK);
phi_star = nan(1, nK);

%% ===== Main loop: compute equilibrium and J_plus for each K =====
for idx = 1:nK
    K = K_vals(idx);
    params = params_template; 
    params.K = K;

    % Vector field in Theta32 frame:
    % x = [R31; R32; phi],  returns [dR31; dR32; dphi]
    f_ref = @(x) rhs_ref_frame(x, params);

    % ---- Primary attempt: fsolve; fallback: simple Newton (no angle wrapping) ----
    x_star = x_guess;
    success = false;
    try
        opts = optimoptions('fsolve','Display','off','FunctionTolerance',tol_root,...
                            'StepTolerance',1e-12,'OptimalityTolerance',1e-12);
        x_star = fsolve(@(y) f_ref(y), x_guess, opts);
        success = true;
    catch
        % Newton fallback with numerical Jacobian
        x = x_guess;
        for it = 1:max_iter_nr
            F = f_ref(x);
            if norm(F) < tol_root, break; end

            Jnum = jacobian_fd_3d(f_ref, x, [1e-6; 1e-6; 1e-6]);

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
        warning('Failed to find a reliable equilibrium at K=%.3f. Skipping.', K);
        % Keep the previous successful x_guess for next continuation step
        continue
    end

    % Record equilibrium
    R31s = x_star(1); 
    R32s = x_star(2); 
    phis = x_star(3);

    R31_star(idx) = R31s; 
    R32_star(idx) = R32s; 
    phi_star(idx) = phis;

    % ---- Compute the 2×2 Jacobian of z_{3+} subsystem (fixing R32 = R32*) ----
    g_plus = @(u) sub_rhs_plus(u, R32s, params);    % u=[R31; phi] -> [dR31; dphi]
    Jplus  = jacobian_fd_2d(g_plus, [R31s; phis], [fd_h_R; fd_h_phi]);
    lams   = eig(Jplus);

    % Save eigenvalues for plotting
    eig1(idx) = lams(1);
    eig2(idx) = lams(2);
    maxReal(idx) = max(real(lams));

    % Continuation: use current equilibrium as initial guess for next K
    x_guess = x_star;
end

%% ===== Plot: Real parts of eigenvalues vs K =====
figure('Color','w'); hold on; grid on; box on;
plot(K_vals, real(eig1), '-', 'LineWidth', 1.8);
plot(K_vals, real(eig2), '-', 'LineWidth', 1.8);
plot(K_vals, maxReal, '--', 'LineWidth', 1.5);   % Track most unstable direction
yline(0,'k:','LineWidth',1.0);

xlabel('K'); 
ylabel('Real part of eigenvalues');
title('Real(\lambda) of J_{plus} versus K  (Theta_{32} frame, R_{32} fixed)');
legend({'Re(\lambda_1)','Re(\lambda_2)','max Re(\lambda)'}, 'Location','best');

% Identify approximate critical K where max Re(lambda) ≈ 0
[~, iCrit] = min(abs(maxReal));
if isfinite(K_vals(iCrit))
    xline(K_vals(iCrit),'r:','LineWidth',1.2);
    text(K_vals(iCrit), 0, sprintf('  K_c \\approx %.3f', K_vals(iCrit)), ...
        'Color',[0.8 0 0], 'VerticalAlignment','bottom');
end

%% ===== Plot equilibrium coordinates for inspection =====
figure('Color','w');
tiledlayout(3,1,'Padding','compact','TileSpacing','compact');

nexttile; 
plot(K_vals, R31_star,'LineWidth',1.6); 
grid on; 
ylabel('R_{31}^*'); 
title('Equilibria vs K');

nexttile; 
plot(K_vals, R32_star,'LineWidth',1.6); 
grid on; 
ylabel('R_{32}^*');

nexttile; 
plot(K_vals, phi_star,'LineWidth',1.6); 
grid on; 
ylabel('\phi^*'); 
xlabel('K');

%% ===== Vector field & numerical Jacobians (no angle wrapping) =====
function F = rhs_ref_frame(x, params)
% Vector field in the Theta32 frame.
% State: x = [R31; R32; phi], where Theta32=0, Theta31=phi.

    R31 = x(1); 
    R32 = x(2); 
    phi = x(3);

    th31 = phi; 
    th32 = 0;

    qfun = params.qfun;
    g = params.gamma;

    C = @(th) g(1)*exp(1i*(th/3)) ...
            + g(2)*exp(1i*(th/3 + 2*pi/3)) ...
            + g(3)*exp(1i*(th/3 + 4*pi/3));

    z11 = qfun(R31) * C(th31);
    z12 = qfun(R32) * C(th32);

    z1  = 0.5*(z11 + z12);
    R1  = abs(z1); 
    th1 = angle(z1);

    K = params.K; 
    Delta = params.Delta; 
    omega0 = params.omega0;

    % ODEs
    dR31   =  3*K/2 * R1^3 * cos(3*th1 - th31) * (1 - R31^2) - 3*R31*Delta;
    dth31  =  3*K/2 * R1^3 * sin(3*th1 - th31) * (R31 + 1/R31) - 3*omega0;

    dR32   =  3*K/2 * R1^3 * cos(3*th1 - th32) * (1 - R32^2) - 3*R32*Delta;
    dth32  =  3*K/2 * R1^3 * sin(3*th1 - th32) * (R32 + 1/R32) + 3*omega0;

    dphi   = dth31 - dth32;   % Relative phase dynamics

    F = [dR31; dR32; dphi];
end

function G = sub_rhs_plus(u, R32fix, params)
% Subsystem for z_{3+}: variables u=[R31; phi], with R32 fixed at R32fix.
    R31 = u(1); 
    phi = u(2);
    x = [R31; R32fix; phi];
    F = rhs_ref_frame(x, params);
    G = [F(1); F(3)];   % Return only [dR31; dphi]
end

function J = jacobian_fd_2d(f, u, h)
% Centered finite-difference Jacobian (2D)
    if numel(h)==1, h = h*ones(2,1); end
    J = zeros(2,2);
    for j=1:2
        up = u; um = u;
        up(j)=up(j)+h(j); 
        um(j)=um(j)-h(j);
        fp = f(up); 
        fm = f(um);
        J(:,j) = (fp - fm)/(2*h(j));
    end
end

function J = jacobian_fd_3d(f, x, h)
% Centered finite-difference Jacobian (3D)
    if numel(h)==1, h = h*ones(3,1); end
    J = zeros(3,3);
    for j=1:3
        xp = x; xm = x;
        xp(j)=xp(j)+h(j); 
        xm(j)=xm(j)-h(j);
        fp = f(xp); 
        fm = f(xm);
        J(:,j) = (fp - fm)/(2*h(j));
    end
end
