% realtime_euler_scan_style_anim_refTheta32_with_equilibrium_and_Jplus.m
% This script performs:
%   (1) Explicit Euler time integration with real-time animation in the
%       rotating frame defined by Theta32.
%   (2) After the animation ends, it searches for an equilibrium
%       (R31*, R32*, phi*) in the Theta32 frame.
%   (3) The detected equilibrium is marked in the projection plot.
%   (4) A 2×2 local Jacobian of z_{3+} (with respect to [R31, phi], holding
%       R32 = R32*) is computed and its eigenvalues are reported.

% clear; clc;

%% ===== Parameters =====
gamma1 = 0.9;  gamma2 = 0.1;  gamma3 = 0;
Delta  = 1.0;
dt     = 1e-5;
Tmax   = 100;
K      = 190;
omega0 = 70;

plotEvery = 10;

%% ===== Load q(|z|) Data (Spline Interpolation with Extrapolation) =====
load('R1_910.mat');   R1_data = solutions(:);
load('R3_910.mat');   R3_data = solutions3(:);
qfun = @(r) interp1(R3_data, R1_data, r, 'spline');   % No [0,1] clipping

%% ===== Definition of C(theta) =====
C = @(th) gamma1*exp(1i*(th/3)) ...
        + gamma2*exp(1i*(th/3 + 2*pi/3)) ...
        + gamma3*exp(1i*(th/3 + 4*pi/3));

%% ===== Initial State =====
% R31 = 1;
Theta31 = 0;
% R32 = 1;
Theta32 = 3.14*0.9;

% Compute z_1, R_1, Theta_1 at the initial state
z11 = qfun(R31) * C(Theta31);
z12 = qfun(R32) * C(Theta32);
z1  = 0.5*(z11 + z12);
R1  = abs(z1);
Theta1 = angle(z1);

%% ===== Set up Figure =====
figure('Color','w'); hold on;
th = linspace(0,2*pi,400);
plot(cos(th), sin(th), 'k-', 'LineWidth',1.2);     % Unit circle
axis equal; axis([-1 1 -1 1]); grid on; box on;
xlabel('Re(z)'); ylabel('Im(z)');
title(sprintf(['Explicit Euler (2nd dynamics)  [frame: \\Theta_{32}]  ', ...
               'K=%.2f, \\omega_0=%.2f, \\Delta=%.1f, dt=%.1e'], ...
               K, omega0, Delta, dt));

traj_plus  = animatedline('Color',[0.85 0.1 0.1],'LineWidth',1.6);
traj_minus = animatedline('Color',[0.1 0.1 0.9],'LineWidth',1.6);

pt_plus  = plot(NaN,NaN,'ro','MarkerFaceColor',[1 0.6 0.6]);
pt_minus = plot(NaN,NaN,'bo','MarkerFaceColor',[0.6 0.6 1]);

legend({'|z|=1','z_{3+} in frame of \Theta_{32}','z_{3-} in frame of \Theta_{32}'}, ...
        'Location','southoutside');

%% ===== Explicit Euler Integration with Real-Time Animation =====
t = 0; k = 0;
while t < Tmax
    k = k + 1;

    % ---- Compute z_1, R_1, Theta_1 for current state ----
    z11 = qfun(R31) * C(Theta31);
    z12 = qfun(R32) * C(Theta32);
    z1  = 0.5*(z11 + z12);
    R1  = abs(z1);
    Theta1 = angle(z1);

    % ---- Second-order Kuramoto-type dynamics ----
    % No R→0 protection and no projection of R onto [0,1]
    dR31     =  3*K/2 * R1^3 * cos(3*Theta1 - Theta31) * (1 - R31^2) - 3*R31*Delta;
    dTheta31 =  3*K/2 * R1^3 * sin(3*Theta1 - Theta31) * (R31 + 1/R31) - 3*omega0;
    dR32     =  3*K/2 * R1^3 * cos(3*Theta1 - Theta32) * (1 - R32^2) - 3*R32*Delta;
    dTheta32 =  3*K/2 * R1^3 * sin(3*Theta1 - Theta32) * (R32 + 1/R32) + 3*omega0;

    % ---- Explicit Euler step ----
    R31     = R31     + dt*dR31;
    Theta31 = Theta31 + dt*dTheta31;
    R32     = R32     + dt*dR32;
    Theta32 = Theta32 + dt*dTheta32;
    t = t + dt;

    % ---- Projection in Theta32 frame: subtract Theta32 ----
    th1_rel = (Theta31 - Theta32)/3;  % Relative angle of z_{3+}
    th2_rel = 0;                      % z_{3-} sits on the real axis

    x1 = R31*cos(th1_rel);  y1 = R31*sin(th1_rel);
    x2 = R32*cos(th2_rel);  y2 = R32*sin(th2_rel);  % equals (R32,0)

    if mod(k, plotEvery) == 0
        addpoints(traj_plus,  x1, y1);
        addpoints(traj_minus, x2, y2);
        set(pt_plus,  'XData', x1, 'YData', y1);
        set(pt_minus, 'XData', x2, 'YData', y2);
        drawnow limitrate;
    end

    % ---- Halt if numerical divergence occurs ----
    if ~isfinite(R31) || ~isfinite(R32) || ~isfinite(Theta31) || ~isfinite(Theta32)
        warning('Numerical divergence/instability detected. Stopping.'); 
        break;
    end
end

%% ===== Search for an Equilibrium in Theta32 Frame =====
% State variable: x = [R31; R32; phi], where phi = Theta31 - Theta32
phi_end = Theta31 - Theta32;
x_end  = [R31; R32; phi_end];

params.gamma = [gamma1, gamma2, gamma3];
params.K = K; params.Delta = Delta; params.omega0 = omega0; params.qfun = qfun;

f_new = @(x) rhs_ref_frame(x, params);   % Vector field in the Theta32 frame

% Equilibrium search using fsolve (fallback: simplified Newton)
x_star = x_end;
try
    opts = optimoptions('fsolve','Display','off','FunctionTolerance',1e-12,...
                        'StepTolerance',1e-12,'OptimalityTolerance',1e-12);
    x_star = fsolve(@(y) f_new(wrap_state(y)), x_end, opts);
    x_star = wrap_state(x_star);
catch
    % Simplified Newton iteration
    for it = 1:30
        F = f_new(x_star);
        if norm(F) < 1e-10, break; end
        J3 = jacobian_fd_3d(f_new, x_star, [1e-6;1e-6;1e-6]);
        dx = - J3 \ F;
        x_star = wrap_state(x_star + dx);
    end
end

R31s = x_star(1); R32s = x_star(2); phis = x_star(3);

%% ===== Plot Detected Equilibrium in the Theta32 Frame =====
x1s = R31s*cos(phis);  y1s = R31s*sin(phis);     % z_{3+}*
x2s = R32s;             y2s = 0;                 % z_{3-}* on x-axis

plot(x1s, y1s, 'r*', 'MarkerSize', 10, 'LineWidth', 1.6);
plot(x2s, y2s, 'b*', 'MarkerSize', 10, 'LineWidth', 1.6);
text(x1s, y1s, '  z_{3+}^*', 'Color',[0.8 0 0],   'FontWeight','bold');
text(x2s, y2s, '  z_{3-}^*', 'Color',[0   0 0.8], 'FontWeight','bold');
drawnow;

%% ===== Compute Local 2×2 Jacobian for z_{3+} Subsystem =====
% Variables: u = [R31; phi], with R32 fixed at R32*
g_plus = @(u) sub_rhs_plus(u, R32s, params);   % Returns [dR31; dphi]
Jplus  = jacobian_fd_2d(g_plus, [R31s; phis], [1e-6; 1e-6]);
eigJp  = eig(Jplus);

fprintf('\n[Theta32 frame] Equilibrium found: R31*=%.6g, R32*=%.6g, phi*=%.6g (rad)\n',...
        R31s, R32s, phis);
fprintf('Local 2x2 Jacobian J_plus = \n'); disp(Jplus);
fprintf('eig(J_plus) = [%.6g%+.6gi, %.6g%+.6gi]\n', ...
        real(eigJp(1)), imag(eigJp(1)), real(eigJp(2)), imag(eigJp(2)));

%% ===== Auxiliary Functions =====
function F = rhs_ref_frame(x, params)
% Vector field in Theta32 frame.
% State: x = [R31; R32; phi], with Theta32 = 0 and Theta31 = phi.

    R31 = x(1); R32 = x(2); phi = x(3);
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

    K = params.K; Delta = params.Delta; omega0 = params.omega0;

    dR31   =  3*K/2 * R1^3 * cos(3*th1 - th31) * (1 - R31^2) - 3*R31*Delta;
    dth31  =  3*K/2 * R1^3 * sin(3*th1 - th31) * (R31 + 1/R31) - 3*omega0;
    dR32   =  3*K/2 * R1^3 * cos(3*th1 - th32) * (1 - R32^2) - 3*R32*Delta;
    dth32  =  3*K/2 * R1^3 * sin(3*th1 - th32) * (R32 + 1/R32) + 3*omega0;

    dphi   = dth31 - dth32;    % Relative phase dynamics

    F = [dR31; dR32; dphi];
end

function G = sub_rhs_plus(u, R32fix, params)
% z_{3+} subsystem in Theta32 frame, with u = [R31; phi].
% R32 is held fixed at R32fix.

    R31 = u(1); 
    phi = u(2);
    x = [R31; R32fix; phi];
    F = rhs_ref_frame(x, params);
    G = [F(1); F(3)];   % Return [dR31; dphi]
end

function xw = wrap_state(x)
% Wrap state vector (currently only phi wrapping possible if needed).
    xw = x;
    xw(3) = xw(3);   % phi unchanged here (placeholder for future wrapping)
end

function J = jacobian_fd_3d(f, x, h)
% 3D finite-difference Jacobian.
    if numel(h)==1, h = h*ones(3,1); end
    J = zeros(3,3);
    for j=1:3
        xp = x; xm = x;
        xp(j)=xp(j)+h(j); 
        xm(j)=xm(j)-h(j);
        xp = wrap_state(xp); 
        xm = wrap_state(xm);
        fp = f(xp); 
        fm = f(xm);
        J(:,j) = (fp - fm)/(2*h(j));
    end
end

function J = jacobian_fd_2d(f, u, h)
% 2D finite-difference Jacobian.
    if numel(h)==1, h = h*ones(2,1); end
    J = zeros(2,2);
    for j=1:2
        up = u; um = u;
        up(j)=up(j)+h(j); 
        um(j)=um(j)-h(j);
        up(2)=up(2); 
        um(2)=um(2);   % phi wrapping placeholder
        fp = f(up); 
        fm = f(um);
        J(:,j) = (fp - fm)/(2*h(j));
    end
end
