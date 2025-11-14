% ================== Order-parameter angular velocity Ω(K) ==================
% Model: generalized Kuramoto phase oscillators with m-th harmonic coupling
%
%     dθ_i/dt = ω_i + K * R^m * sin[m(Θ − θ_i)]
%
% where the complex order parameter is
%
%     Z = R e^{iΘ} = (1/N) ∑_j e^{iθ_j}.
%
% For m = 1,2,3,4 the script simulates the dynamics over a grid of K-values
% and extracts the collective angular velocity Ω = dΘ/dt from the long-time
% tail of the trajectory.  Values with tail-averaged |Z| < Rmin are omitted.
% ===========================================================================

clear; clc; close all;
rng(1);

% ----------------------------- Global parameters ---------------------------
N      = 1000;        % Number of oscillators
sigma  = 1.0;         % Lorentzian width for ω_i
T      = 20.0;        % Total simulation time
dt     = 1e-3;        % Time step
nSteps = round(T/dt);
tvec   = (0:nSteps-1)*dt;

Rmin   = 0.1;         % Minimum |Z| to accept Ω(K) data

% K-ranges for each harmonic order m
Kgrid.m1 = linspace(0, 30, 30);
Kgrid.m2 = linspace(0, 30, 30);
Kgrid.m3 = linspace(0, 30, 30);
Kgrid.m4 = linspace(0, 30, 30);

mList = [1 2 3 4];

% Initial phase distribution: populations at 0, 2π/3, 4π/3 with ratio 9:1:0
ratios = [0.9 0.1 0.0];
pts    = [0, 2*pi/3, 4*pi/3];

% Tail window used for steady-state statistics
tailFrac = 0.3;
tailIdx  = round((1-tailFrac)*nSteps):nSteps;

% ---------------------- Natural frequencies (Lorentzian) -------------------
U = rand(N,1);
omega = sigma * tan(pi*(U - 0.5));

% ---------------------------- Initial phases -------------------------------
theta0 = zeros(N,1);
counts = round(N*ratios);
counts(end) = N - sum(counts(1:end-1));    % Enforce total = N
idx = 1;
for g = 1:3
    if counts(g) > 0
        theta0(idx:idx+counts(g)-1) = pts(g);
        idx = idx + counts(g);
    end
end

% ======================= Main loop: m = 1,2,3,4 ============================
figure('Color','w'); hold on; grid on; box on;
colors  = lines(numel(mList));
markers = {'o','s','^','d'};
h = gobjects(1,numel(mList));
legendNames = cell(1,numel(mList));

for im = 1:numel(mList)
    m = mList(im);

    % Select K-range for this m
    switch m
        case 1, Kvals = Kgrid.m1;
        case 2, Kvals = Kgrid.m2;
        case 3, Kvals = Kgrid.m3;
        case 4, Kvals = Kgrid.m4;
    end

    Omega_vs_K = nan(size(Kvals));
    Rtail_mean = nan(size(Kvals));

    fprintf('Simulating m = %d ...\n', m);

    for ik = 1:numel(Kvals)
        K = Kvals(ik);
        theta = theta0;

        Ztail = zeros(1, numel(tailIdx));
        tcnt  = 0;

        for it = 1:nSteps

            % Current order parameter Z
            Z = mean(exp(1i*theta));
            R = abs(Z);  Th = angle(Z);

            % Phase velocity
            f = @(th) omega + K*(R^m).*sin(m*(Th - th));

            % RK4 step
            k1 = f(theta);

            Z2 = mean(exp(1i*(theta + 0.5*dt*k1)));
            f2 = @(th) omega + K*(abs(Z2)^m).*sin(m*(angle(Z2) - th));
            k2 = f2(theta + 0.5*dt*k1);

            Z3 = mean(exp(1i*(theta + 0.5*dt*k2)));
            f3 = @(th) omega + K*(abs(Z3)^m).*sin(m*(angle(Z3) - th));
            k3 = f3(theta + 0.5*dt*k2);

            Z4 = mean(exp(1i*(theta + dt*k3)));
            f4 = @(th) omega + K*(abs(Z4)^m).*sin(m*(angle(Z4) - th));
            k4 = f4(theta + dt*k3);

            theta = theta + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
            theta = angle(exp(1i*theta));   % Wrap to (-π, π]

            % Record tail window only
            if it >= tailIdx(1)
                tcnt = tcnt + 1;
                Ztail(tcnt) = mean(exp(1i*theta));
            end
        end

        % Tail statistics
        Rt = abs(Ztail);
        Rtail_mean(ik) = mean(Rt);

        if Rtail_mean(ik) >= Rmin
            ph = unwrap(angle(Ztail));
            tt = tvec(tailIdx);
            p  = polyfit(tt, ph, 1);
            Omega_vs_K(ik) = p(1);
        else
            Omega_vs_K(ik) = NaN;
        end
    end

    % Plot Ω(K) for this m
    h(im) = plot(Kvals, Omega_vs_K, ...
        'LineWidth', 2, 'MarkerSize', 5, ...
        'Color', colors(im,:), 'Marker', markers{im});
    legendNames{im} = sprintf('m = %d', m);
end

% ------------------------------ Final figure -------------------------------
yline(0,'k--','LineWidth',1);
xlabel('K');
ylabel('\Omega = d/dt arg(Z)');
title(sprintf('\\Omega(K) for m = 1,2,3,4  (valid only for \\langle|Z|\\rangle_{tail} \\ge %.2f)', Rmin));
legend(h, legendNames, 'Location','best');
hold off;
