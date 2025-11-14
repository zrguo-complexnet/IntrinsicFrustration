% Parallel sweep over K: from 110 down to 95. For each K, compute the
% corresponding Floquet multiplier. Failed K values are automatically skipped.

clc; clear;

%% ===== Fixed parameters (consistent with the main script) =====
gamma1 = 0.9; %#ok<NASGU>
gamma2 = 0.1; %#ok<NASGU>
gamma3 = 0.0; %#ok<NASGU>

omega0 = 100;
Delta  = 1;
dt     = 1e-5;
Tmax   = 5;

% Initial conditions (same for every K)
R31_0     = 1.0;  
Theta31_0 = 0.0;
R32_0     = 1.0;  
Theta32_0 = 0.0;

% R1(R3) data and interpolant (griddedInterpolant/interp1: parallel-safe)
load('R1_910.mat');  % solutions -> R1_data
load('R3_910.mat');  % solutions3 -> R3_data
R1_data = solutions;
R3_data = solutions3;

R1_from_R3 = @(R3) interp1(R3_data, R1_data, R3, 'spline', 'extrap');

% Poincaré section & Floquet parameters
section_dir = -1;       % direction of crossing
M_settle    = 6;        % number of consecutive section hits required for convergence
tol_section = 1e-3;     % tolerance for convergence on the Poincaré section
perturb_var = 'R31';    % variable used for perturbation
eps_perturb = 1e-1;     % perturbation magnitude
verbose     = false;    % disable printing in parallel mode

% K sweep (descending order)
K_vec = 110:-0.2:95;
nK = numel(K_vec);

% Containers for results (parallel-friendly)
mu_all = nan(nK,1);
K_ok   = nan(nK,1);

%% ===== Parallel pool =====
if isempty(gcp('nocreate'))
    parpool;    % start pool with maximum available workers
end

%% ===== Parallel loop =====
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
        % If anything fails at this K, simply skip it.
    end
end

%% ===== Collect successful results and plot =====
mask = isfinite(mu_all);
K_keep  = K_ok(mask);
mu_keep = mu_all(mask);

figure('Color','w');
plot(K_keep, mu_keep, 'o-', ...
     'LineWidth',1.5, 'MarkerSize',5);
set(gca,'XDir','reverse');   % display K from 110 down to 95
xlabel('K'); 
ylabel('\mu (Floquet multiplier)');
title('Floquet Multiplier vs. K (parallel sweep)');
grid on;


%% ========================================================================
%% ========================= Helper Functions =============================
%% ========================================================================

function [ok, mu] = run_one_K(K, omega0, Delta, dt, Tmax, ...
                              R31, Theta31, R32, Theta32, ...
                              R1_from_R3, section_dir, M_settle, tol_section, ...
                              perturb_var, eps_perturb, verbose)
% Full workflow for a single K:
%   1) integrate until convergence onto a stable limit cycle,
%   2) detect Poincaré section crossings,
%   3) check stability across consecutive section hits,
%   4) compute the Floquet multiplier via central-difference perturbation.

    ok = false; 
    mu = NaN;

    Nmax = round(Tmax/dt);
    phi_prev = wrapToPi(Theta31 - Theta32);
    t = 0;

    R_section = [];
    state_at_section = [];
    t_at_section = [];

    % ---- Main integration loop, including Poincaré section detection ----
    for k = 1:Nmax
        t = t + dt;

        % RK4 time-stepping
        [k1_R31, k1_T31, k1_R32, k1_T32] = rhs(R31,Theta31,R32,Theta32,K,omega0,Delta,R1_from_R3);
        [k2_R31, k2_T31, k2_R32, k2_T32] = rhs(R31+0.5*dt*k1_R31, Theta31+0.5*dt*k1_T31, ...
                                              R32+0.5*dt*k1_R32, Theta32+0.5*dt*k1_T32, ...
                                              K,omega0,Delta,R1_from_R3);
        [k3_R31, k3_T31, k3_R32, k3_T32] = rhs(R31+0.5*dt*k2_R31, Theta31+0.5*dt*k2_T31, ...
                                              R32+0.5*dt*k2_R32, Theta32+0.5*dt*k2_T32, ...
                                              K,omega0,Delta,R1_from_R3);
        [k4_R31, k4_T31, k4_R32, k4_T32] = rhs(R31+dt*k3_R31, Theta31+dt*k3_T31, ...
                                              R32+dt*k3_R32, Theta32+dt*k3_T32, ...
                                              K,omega0,Delta,R1_from_R3);

        % Update state
        R31     = R31     + dt/6*(k1_R31 + 2*k2_R31 + 2*k3_R31 + k4_R31);
        Theta31 = Theta31 + dt/6*(k1_T31 + 2*k2_T31 + 2*k3_T31 + k4_T31);
        R32     = R32     + dt/6*(k1_R32 + 2*k2_R32 + 2*k3_R32 + k4_R32);
        Theta32 = Theta32 + dt/6*(k1_T32 + 2*k2_T32 + 2*k3_T32 + k4_T32);

        % Poincaré section detection
        phi = wrapToPi(Theta31 - Theta32);
        crossed = (section_dir > 0 && (phi_prev < 0 && phi >= 0 && cos(phi/3) > 0)) || ...
                  (section_dir < 0 && (phi_prev > 0 && phi <= 0 && cos(phi/3) > 0));

        if crossed
            % Linear interpolation to estimate exact section crossing
            alpha = abs(phi_prev) / (abs(phi_prev) + abs(phi) + eps);

            R31_sec     = R31     - alpha * (R31     - (R31 - dt/6*(k1_R31 + 2*k2_R31 + 2*k3_R31 + k4_R31)));
            Theta31_sec = Theta31 - alpha * (Theta31 - (Theta31 - dt/6*(k1_T31 + 2*k2_T31 + 2*k3_T31 + k4_T31)));
            R32_sec     = R32     - alpha * (R32     - (R32 - dt/6*(k1_R32 + 2*k2_R32 + 2*k3_R32 + k4_R32)));
            Theta32_sec = Theta32 - alpha * (Theta32 - (Theta32 - dt/6*(k1_T32 + 2*k2_T32 + 2*k3_T32 + k4_T32)));
            t_sec       = t - alpha*dt;

            % Select variable to use for Floquet computation
            switch perturb_var
                case 'R31', R_on_sec = R31_sec;
                case 'R32', R_on_sec = R32_sec;
                otherwise
                    error('perturb_var must be ''R31'' or ''R32''.');
            end

            % Store section-crossing information
            R_section(end+1,1) = R_on_sec; %#ok<AGROW>
            state_at_section(end+1,:) = [R31_sec,Theta31_sec,R32_sec,Theta32_sec]; %#ok<AGROW>
            t_at_section(end+1,1) = t_sec; %#ok<AGROW>

            % Check if section hits have converged → stable limit cycle detected
            if numel(R_section) >= M_settle+1
                recent = R_section(end-M_settle:end);
                if max(abs(diff(recent))) < tol_section

                    % Extract limit cycle section state
                    R_star = R_section(end);
                    x_star = state_at_section(end,:);

                    % Compute Floquet multiplier via central difference
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


%% ========================= System RHS (no wrapping changes) =========================
function [dR31,dT31,dR32,dT32] = rhs(R31,Theta31,R32,Theta32,K,omega0,Delta,R1_from_R3)
% Standard ODE right-hand side for the (R31,Theta31,R32,Theta32) system.

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
    R1 = abs(z1);
    Theta1 = angle(z1);

    R31 = max(R31, 1e-12);
    R32 = max(R32, 1e-12);

    dR31 = (3*K/2) * R1^3 * cos(3*Theta1 - Theta31) * (1 - R31^2) - 3*Delta*R31;
    dT31 = (3*K/2) * R1^3 * sin(3*Theta1 - Theta31) * (R31 + 1/R31) - 3*omega0;

    dR32 = (3*K/2) * R1^3 * cos(3*Theta1 - Theta32) * (1 - R32^2) - 3*Delta*R32;
    dT32 = (3*K/2) * R1^3 * sin(3*Theta1 - Theta32) * (R32 + 1/R32) + 3*omega0;
end


%% ========================= Floquet Multiplier (central difference) =========================
function mu = computeFloquet(x_star, perturb_var, eps0, ...
                             K,omega0,Delta,dt,section_dir,R1_from_R3, ...
                             ~,~,~,~,~,verbose)
% Compute a scalar Floquet multiplier associated with perturbations along
% the chosen direction (R31 or R32), using a central finite-difference
% approximation after one section-to-section return.

    eps = eps0;

    % Perturbed initial points
    x_plus  = x_star; 
    x_minus = x_star;

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

    % Return map evaluations
    Rp_next = advance_to_next_section(x_plus(1),x_plus(2),x_plus(3),x_plus(4), ...
        K,omega0,Delta,dt,section_dir,R1_from_R3);

    Rm_next = advance_to_next_section(x_minus(1),x_minus(2),x_minus(3),x_minus(4), ...
        K,omega0,Delta,dt,section_dir,R1_from_R3);

    % Central-difference Floquet multiplier
    switch perturb_var
        case 'R31', mu = (Rp_next(1) - Rm_next(1)) / (2*eps);
        case 'R32', mu = (Rp_next(3) - Rm_next(3)) / (2*eps);
    end

    if verbose
        fprintf('  (central difference collected)\n');
    end
end


%% ========================= Section-to-Section Map =========================
function x_next = advance_to_next_section(R31,Theta31,R32,Theta32, ...
        K,omega0,Delta,dt,section_dir,R1_from_R3)
% Integrate forward until the next crossing of the Poincaré section
% (same definition as in the main loop). Linear interpolation ensures
% sub-step accuracy for the section return.

    % First, take one step off the section
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

    % Continue integration until section crossing occurs
    while true
        % RK4 loop
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
            % Linear interpolation for sub-step accuracy
            alpha = abs(phi_prev) / (abs(phi_prev) + abs(phi) + eps);

            R31_sec     = R31_new     - alpha * (R31_new     - R31);
            Theta31_sec = Theta31_new - alpha * (Theta31_new - Theta31);
            R32_sec     = R32_new     - alpha * (R32_new     - R32);
            Theta32_sec = Theta32_new - alpha * (Theta32_new - Theta32);

            x_next = [R31_sec,Theta31_sec,R32_sec,Theta32_sec];
            return;
        end

        % Update state
        R31=R31_new; 
        Theta31=Theta31_new; 
        R32=R32_new; 
        Theta32=Theta32_new; 
        phi_prev=phi;
    end
end

function a = wrapToPi(a)
% Wrap angle into the symmetric interval of width 2π.
    a = mod(a + 3*pi, 6*pi) - 3*pi;
end
