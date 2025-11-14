% HC_DrawManifold.m
% Purpose: Numerically detect saddle-type equilibria for a given coupling strength K,
%          and classify their stability types.
% Dependencies: R1_910.mat, R3_910.mat (consistent with the user’s existing code)

% clear; clc; close all;

%% ===== Parameter Settings =====
gamma1 = 0.9;  gamma2 = 0.1;  gamma3 = 0.0;
Delta  = 1.0;
omega0 = 70;

% Select the value of K to be examined.
% It is recommended to test values on both sides of the saddle-node bifurcation threshold.
K = 177;

% Numerical parameters for root-finding / finite differences
tol_root    = 1e-10;            % Termination threshold: ||F||_2
max_newton  = 40;               % Maximum Newton iterations
fd_step_3d  = [1e-6;1e-6;1e-6]; % Finite-difference step sizes for Jacobian
alpha_list  = [1e-3, 5e-3, 1e-2, 5e-2]; % Trial step sizes along eigenvectors
n_multistart= 80;               % Fallback multi-start sample size

% Initial condition (first find a stable equilibrium)
x_init = [0.9175; 0.9835; -18.684460444436220-2.286780062262495];    % [R31; R32; phi] in Theta32 frame, no wrapping

%% ===== Data and Vector Field =====
load('R1_910.mat');   R1_data = solutions(:);
load('R3_910.mat');   R3_data = solutions3(:);

% Cubic self-consistency map q(r), evaluated by spline interpolation.
qfun = @(r) interp1(R3_data, R1_data, r, 'spline');   % Extrapolation allowed

params.gamma = [gamma1,gamma2,gamma3];
params.K = K; params.Delta = Delta; params.omega0 = omega0; params.qfun = qfun;

F    = @(x) rhs_ref_frame(x, params);     % Vector field F(x) = [dR31; dR32; dphi]
Jnum = @(x) jacobian_fd_3d(F, x, fd_step_3d);

%% ===== Step 1: Find a Stable Equilibrium (Anchor Point) =====
[x_stb, ok_stb] = solve_equilibrium(F, Jnum, x_init, tol_root, max_newton);

% Use fsolve as backup if Newton fails
if ~ok_stb
    try
        opts = optimoptions('fsolve','Display','off','FunctionTolerance',tol_root,...
                            'StepTolerance',1e-12,'OptimalityTolerance',1e-12);
        x_stb = fsolve(@(y) F(y), x_init, opts);
        ok_stb = true;
    catch
        ok_stb = false;
    end
end
assert(ok_stb, 'Failed to converge to a stable equilibrium. Adjust x_init or K.');

% Stability classification
J_stb = Jnum(x_stb);
lam_stb = eig(J_stb);

fprintf('Stable equilibrium x_stb = [%.6g, %.6g, %.6g]^T\n', x_stb(1), x_stb(2), x_stb(3));
fprintf('eig(J_stb) = %s\n', sprintf('[%.3e%+.3ei] ', real(lam_stb(1)), imag(lam_stb(1))));
stb_type = classify_eig(lam_stb);
fprintf('Type (stable/saddle/source/critical): %s\n\n', stb_type);

%% ===== Step 2: Search for Saddle Points by Shooting Along the Smallest-Eigenvalue Direction =====
[~, idx_min] = min(abs(real(lam_stb)));   % Select eigenvalue closest to 0 in real part
[v_min, ~]   = eigvec_of(J_stb, idx_min);

cands = [];  % Candidate equilibria
for sgn = [-1, 1]
    for a = alpha_list
        x0 = x_stb + sgn * a * real(v_min);
        [x_eq, ok] = solve_equilibrium(F, Jnum, x0, tol_root, max_newton);
        if ok
            cands = [cands, x_eq]; %#ok<AGROW>
        end
    end
end

% Merge approximately identical solutions
eqs = unique_points(cands, 1e-6);

%% ===== Fallback: Multi-Start Random Search if Vector-Shooting Finds No New Points =====
if isempty(eqs)
    rng(0);
    cands2 = [];

    % Randomized sampling box (adjust according to system scale)
    R31_box = [0.1, 1.5];
    R32_box = [0.1, 1.5];
    phi_box = [-20, 20];   % No wrapping

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

%% ===== Classify All Found Equilibria and Extract Saddle Points =====
saddles = [];
fprintf('=== List of Detected Equilibria (with Stability Types) ===\n');

if isempty(eqs)
    fprintf('No additional equilibria detected.\n');
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
    fprintf('\nNo saddle points found. You may try adjusting K near criticality, tuning alpha_list, or using denser multi-start sampling.\n');
else
    fprintf('\nDetected %d saddle point(s).\n', size(saddles,2));
end

%% ===== (Optional) Simple Phase-Space Visualization in the Theta32 Frame =====
if ~isempty(saddles)
    figure('Color','w'); hold on; axis equal; grid on; box on;
    th = linspace(0,2*pi,400); plot(cos(th),sin(th),'k-'); % Unit circle

    % Stable equilibria
    plot(x_stb(1)*cos(x_stb(3)), x_stb(1)*sin(x_stb(3)), 'go','MarkerFaceColor',[0.6 1 0.6],'DisplayName','stable z_{3+}^*');
    plot(x_stb(2), 0, 'g^','MarkerFaceColor',[0.6 1 0.6],'DisplayName','stable z_{3-}^*');

    % Saddle equilibria (red markers)
    for i=1:size(saddles,2)
        xs = saddles(:,i);
        plot(xs(1)*cos(xs(3)), xs(1)*sin(xs(3)), 'rd','MarkerFaceColor',[1 0.6 0.6],'DisplayName','saddle z_{3+}^*');
        plot(xs(2), 0, 'rs','MarkerFaceColor',[1 0.6 0.6],'DisplayName','saddle z_{3-}^*');
    end

    xlabel('Re(z)'); ylabel('Im(z)');
    title(sprintf('Equilibria in \\Theta_{32} Frame (K=%.3f)', K));
    legend('Location','best');
end

%% ===== Live Euler Integration: Stable/Unstable Manifolds at Saddle Points =====
if ~isempty(saddles)
    % Background figure
    figure('Color','w'); hold on; axis equal; grid on; box on;
    th = linspace(0,2*pi,400); plot(cos(th),sin(th),'k-');
    xlabel('Re(z)'); ylabel('Im(z)');
    title(sprintf('Live Euler Manifolds at Saddles (T=1, dt=1e-4), K=%.3f', K));

    % Reference stable equilibrium
    plot(x_stb(1)*cos(x_stb(3)), x_stb(1)*sin(x_stb(3)), 'go','MarkerFaceColor',[0.6 1 0.6],'DisplayName','stable z_{3+}^*');
    plot(x_stb(2), 0, 'g^','MarkerFaceColor',[0.6 1 0.6],'DisplayName','stable z_{3-}^*');

    % Integration parameters
    dt = 1e-4;   % Step size
    T  = 1.0;    % Total time
    nsteps = round(T/dt);

    eps0   = 1e-4;    % Initial perturbation
    Rmin   = 1e-6;    % Lower bound for R
    Rmax   = 1.0;     % Upper bound for R
    phiMax = 40;      % Upper bound for |phi|

    colU = [0.85 0.1 0.1];  % Unstable (red)
    colS = [0.1 0.2 0.85];  % Stable (blue)

    % Loop through saddle equilibria
    for i = 1:size(saddles,2)
        xs = saddles(:,i);

        % Plot saddle projections
        plot(xs(1)*cos(xs(3)), xs(1)*sin(xs(3)), 'kd','MarkerFaceColor',[1 0.8 0.8],'DisplayName','saddle z_{3+}^*');
        plot(xs(2), 0, 'ks','MarkerFaceColor',[1 0.8 0.8],'DisplayName','saddle z_{3-}^*');

        % Jacobian and eigen-decomposition
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

        % Unstable manifold (Euler forward)
        if ~isempty(v_u)
            for sgn = [-1, 1]
                x0 = xs + sgn*eps0*v_u;
                live_euler_and_plot(F, x0, dt, nsteps, Rmin, Rmax, phiMax, colU);
            end
        end

        % Stable manifold (Euler backward)
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
% Euler integration with continuous visualization.
% Includes projection onto disk (z_{3+}) and real axis (z_{3-}).

    x = x0;

    h1 = animatedline('Color',col,'LineWidth',1.4,'MaximumNumPoints',nsteps,'DisplayName','');
    h2 = animatedline('Color',col,'LineWidth',1.0,'MaximumNumPoints',nsteps,'HandleVisibility','off');

    [xp,yp,xm] = proj_one(x);
    addpoints(h1, xp, yp);
    addpoints(h2, xm, 0);

    for k = 1:nsteps
        fx = F(x);
        x  = x + dt * fx;

        if ~all(isfinite(x)), break; end
        if x(1) < Rmin || x(2) < Rmin, break; end
        if x(1) > Rmax || x(2) > Rmax, break; end
        if abs(x(3)) > phiMax,         break; end

        [xp,yp,xm] = proj_one(x);
        addpoints(h1, xp, yp);
        addpoints(h2, xm, 0);

        drawnow limitrate
    end
end

function [xp,yp,xm] = proj_one(x)
% Projection:
% (R31,R32,phi) → (Re(z_{3+}), Im(z_{3+})) and real-axis position of z_{3-}.
    R31 = x(1); R32 = x(2); phi = x(3);
    xp = R31*cos(phi/3);
    yp = R31*sin(phi/3);
    xm = R32;
end

%% ===== Stable / Unstable Manifold Computation at a Saddle Point =====
function plot_manifolds_at_saddle(xs, F, Jfun, opts, style)
% xs:       Saddle equilibrium [R31; R32; phi]
% F:        Right-hand-side function handle
% Jfun:     Jacobian evaluator
% opts:     Numerical parameters (step sizes, cutoffs, etc.)
% style:    Visualization style

    arguments
        xs (3,1) double
        F function_handle
        Jfun function_handle
        opts.h_forward (1,1) double = 1e-3
        opts.h_backward (1,1) double = -1e-3
        opts.nsteps_u (1,1) double = 8000
        opts.nsteps_s (1,1) double = 4000
        opts.eps0 (1,1) double = 1e-4
        opts.n_directions_s (1,1) double = 24
        opts.Rmin (1,1) double = 1e-6
        opts.Rmax (1,1) double = 3
        opts.phi_max (1,1) double = 40
        opts.keep_every (1,1) double = 5
        style.plot3D (1,1) logical = false
        style.plotProjected (1,1) logical = true
        style.colorU (1,3) double = [0.85 0.1 0.1]
        style.colorS (1,3) double = [0.1 0.2 0.85]
    end

    % Eigen-decomposition
    J = Jfun(xs);
    [V,D] = eig(J);
    lam   = diag(D);

    % Classification
    is_pos = real(lam) >  1e-10;
    is_neg = real(lam) < -1e-10;
    is_zero= ~(is_pos | is_neg);

    Vu = V(:,is_pos);         % Unstable subspace
    Vs = V(:,is_neg | is_zero); % Stable/center subspace

    normalize = @(v) v / max(1e-12, norm(real(v)));
    if ~isempty(Vu),  Vu = real(Vu);  Vu(:,1) = normalize(Vu(:,1)); end
    if ~isempty(Vs),  Vs = real(Vs);  Vs = orth(Vs);               end

    % ===== Unstable manifold (usually 1D) =====
    if ~isempty(Vu)
        vu = Vu(:,1);
        for sgn = [-1, 1]
            x0 = xs + sgn * opts.eps0 * vu;
            Xu = euler_traj(F, x0, opts.h_forward, opts.nsteps_u, opts);
            draw_traj(Xu, style.colorU, style);
        end
    end

    % ===== Stable manifold (1D or 2D) =====
    if size(Vs,2)>=1
        if size(Vs,2)==1
            vs1 = normalize(Vs(:,1));
            for sgn = [-1, 1]
                x0 = xs + sgn * opts.eps0 * vs1;
                Xs = euler_traj(F, x0, opts.h_backward, opts.nsteps_s, opts);
                draw_traj(Xs, style.colorS, style);
            end
        else
            % 2D stable subspace: sample multiple directions on the unit circle
            vs1 = normalize(Vs(:,1));
            vs2 = normalize(Vs(:,2));

            thetas = linspace(0, 2*pi, opts.n_directions_s+1); 
            thetas(end)=[];

            for th = thetas
                dir = cos(th)*vs1 + sin(th)*vs2;
                x0  = xs + opts.eps0 * dir;
                Xs  = euler_traj(F, x0, opts.h_backward, opts.nsteps_s, opts);
                draw_traj(Xs, style.colorS, style);
            end
        end
    end

    % Plot the saddle equilibrium itself
    if style.plotProjected
        [xp,yp,xm] = proj_point(xs);
        plot(xp,yp,'kd','MarkerFaceColor',[1 0.8 0.8]);
        plot(xm,0 ,'ks','MarkerFaceColor',[1 0.8 0.8]);
    end
    if style.plot3D
        plot3(xs(1), xs(2), xs(3), 'kd','MarkerFaceColor',[1 0.8 0.8]);
    end
end

function X = euler_traj(F, x0, h, nsteps, opts)
% Simple forward (h>0) or backward (h<0) Euler integration.

    X = zeros(3, floor(nsteps/opts.keep_every)+1);
    x = x0; kkeep = 1; X(:,kkeep) = x;

    for k = 1:nsteps
        fx = F(x);
        x = x + h * fx;

        if ~all(isfinite(x)), break; end
        if x(1) < opts.Rmin || x(2) < opts.Rmin, break; end
        if x(1) > opts.Rmax || x(2) > opts.Rmax, break; end
        if abs(x(3)) > opts.phi_max, break; end

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
        plot(xm,0*xm,'-','Color',col,'LineWidth',1.2);
    end
    if style.plot3D
        plot3(X(1,:), X(2,:), X(3,:), '-', 'Color', col, 'LineWidth', 1.2);
    end
end

function [xp,yp,xm] = proj_points(X)
% Map (R31,R32,phi) to projections:
% z_{3+} = R31 * exp(i*phi)  → (xp, yp)
% z_{3-} = R32               → xm on real axis
    R31 = X(1,:); R32 = X(2,:); phi = X(3,:);
    xp = R31 .* cos(phi);
    yp = R31 .* sin(phi);
    xm = R32;
end

function [xp,yp,xm] = proj_point(x)
    [xp,yp,xm] = proj_points(x);
end

%% ===== Function Definitions =====
function F = rhs_ref_frame(x, params)
% Vector field in the Theta32 rotating frame.
% State: x = [R31; R32; phi], where Theta32 = 0 and Theta31 = phi.

    R31 = x(1); R32 = x(2); phi = x(3);
    th31 = phi; th32 = 0;

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

    % Dynamics (consistent with user's previous code).
    dR31   =  3*K/2 * R1^3 * cos(3*th1 - th31) * (1 - R31^2) - 3*R31*Delta;
    dth31  =  3*K/2 * R1^3 * sin(3*th1 - th31) * (R31 + 1/R31) - 3*omega0;
    dR32   =  3*K/2 * R1^3 * cos(3*th1 - th32) * (1 - R32^2) - 3*R32*Delta;
    dth32  =  3*K/2 * R1^3 * sin(3*th1 - th32) * (R32 + 1/R32) + 3*omega0;

    dphi   = dth31 - dth32;
    F = [dR31; dR32; dphi];
end

function [x, ok] = solve_equilibrium(F, Jfun, x0, tol_root, maxit)
    x = x0; ok = false;
    for it = 1:maxit
        Fx = F(x);
        if norm(Fx) < tol_root
            ok = true;
            break;
        end
        Jx = Jfun(x);
        dx = - Jx \ Fx;
        x  = x + dx;
        if ~all(isfinite(x)), ok=false; return; end
    end
end

function J = jacobian_fd_3d(f, x, h)
    if numel(h)==1, h = h*ones(3,1); end
    J = zeros(3,3);
    for j = 1:3
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
