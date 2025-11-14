% ==========================================================
% run_highorder_kuramoto_on_flag.m
% High-order Kuramoto dynamics (up to 6-body terms) on a
% precomputed random flag complex.
%
% Frequency assignment:
%   first half:  +omega0 + sigma*rand
%   second half: -omega0 - sigma*rand
%
% r-body coupling:
%   contributions along r-simplices; each node's contribution
%   is normalized by its r-degree (number of incident r-simplices).
%
% Initial phases:
%   symmetric two-group pattern (theta1 / flipped theta1).
% ==========================================================
clear; clc;

%% ========= 1) Select network bundle (.mat) and simulation parameters =========
bundleMat = "";  % If empty: automatically pick latest flag_*/flag_complex_bundle.mat
if bundleMat == "", bundleMat = string(findLatestBundle()); end
if bundleMat == "", error('No flag_* directory with flag_complex_bundle.mat found.'); end
Sdata = load(bundleMat);
fprintf('Loaded: %s\n', bundleMat);

% --- Kuramoto parameters ---
omega0 = 5;        % Base frequency (+/-)
sigma  = 0.1;      % Frequency noise amplitude
Kvec   = zeros(6,1);
% Kvec(2..6) correspond to r-body couplings r = 2..6
Kvec(2) = 0.5;     % 2-body
Kvec(3) = 0;       % 3-body
Kvec(4) = 25;      % 4-body
% Kvec(5) = 0;     % 5-body
% Kvec(6) = 0;     % 6-body

dt = 0.1;          % Time step
T  = 20;           % Total simulation time
rng(2025);         % Seed for frequency noise

%% ========= 2) Extract network data & prepare r-degrees =========
A = Sdata.A;       % n×n adjacency matrix
S = Sdata.S;       % S{r}: list of r-vertex cliques (sorted rows)
f = Sdata.f(:);
N = size(A,1);
Rmax = 4;          % Use interactions up to r ≤ 4 (here 2–4 bodies)

% Node r-degrees: number of r-vertex cliques incident to each node
deg = cell(Rmax,1);
if isfield(Sdata,'deg') && numel(Sdata.deg) >= Rmax
    for r = 2:Rmax
        deg{r} = Sdata.deg{r};
    end
else
    for r = 2:Rmax
        Er = S{r};
        if isempty(Er)
            deg{r} = zeros(N,1);
        else
            deg{r} = accumarray(Er(:), 1, [N,1], @sum, 0);
        end
    end
end

%% ========= 3) Natural frequencies & initial phases =========
halfN  = floor(N/2);
omega  = zeros(N,1);
omega(1:halfN)     =  omega0 + sigma*rand(halfN,1);
omega(halfN+1:end) = -omega0 - sigma*rand(N-halfN,1);

rng(2598);
theta1 = zeros(N/2,1) - 2*pi/3 * (rand(N/2,1) < 0.1);   % First-group initial phases
theta2 = flip(theta1);                                  % Second group: mirrored
theta0 = [theta1; theta2];                              % Full initial condition

%% ========= 4) Time stepping (explicit Euler) =========
nT   = floor(T/dt) + 1;
t    = (0:nT-1).' * dt;
Theta = zeros(N, nT);
Theta(:,1) = theta0;

absz = zeros(1, nT);           % |Z| for the whole population
z_half    = zeros(2, nT);      % Complex order parameters of two halves
absz_half = zeros(2, nT);      % Their moduli

for k = 1:nT-1
    theta  = Theta(:,k);
    dtheta = omega;            % Base frequency term

    % r-body couplings (r = 2..Rmax) along existing cliques,
    % normalized by each node's r-degree
    for r = 2:Rmax
        if Kvec(r) == 0, continue; end
        Er = S{r};
        m  = size(Er,1);
        if m == 0, continue; end

        sum_e    = sum(theta(Er), 2);         % Sum of phases on each r-simplex (m×1)
        idx_all  = Er(:);                     % (m*r)×1 unfolded vertex indices
        sum_rep  = kron(ones(r,1), sum_e);    % Repeat sums for each vertex in simplex
        theta_in = theta(idx_all);            % Local phase at each (node, simplex)

        % Interaction term: sin(∑_{j∈e} θ_j − r θ_i)
        contrib = sin(sum_rep - r * theta_in);

        % Normalize by r-degree of each node
        dr = deg{r}(idx_all);                 % r-degree of each involved node
        contrib = contrib ./ dr;

        % Aggregate contributions per node
        dtheta = dtheta + Kvec(r) * ...
                 accumarray(idx_all, contrib, [N,1], @sum, 0);
    end

    % One explicit Euler step
    Theta(:,k+1) = Theta(:,k) + dt * dtheta;

    % Global order parameter
    z = mean(exp(1i*Theta(:,k+1)));
    absz(1,k+1) = abs(z);

    % Order parameters of the two halves
    idx1 = 1:halfN;
    idx2 = (halfN+1):N;
    z1 = mean(exp(1i*Theta(idx1, k+1)));
    z2 = mean(exp(1i*Theta(idx2, k+1)));

    z_half(1, k+1)    = z1;
    z_half(2, k+1)    = z2;
    absz_half(1, k+1) = abs(z1);
    absz_half(2, k+1) = abs(z2);
end

%% ========= 5) Basic visualization of |z| for the two halves =========
figure('Name','Half-group order parameters');
hold on; grid on;

col1 = [0 0.4470 0.7410];
col2 = [0.8500 0.3250 0.0980];

p1 = plot(t, absz_half(1,:), 'LineWidth', 1.8, 'Color', col1);
p2 = plot(t, absz_half(2,:), 'LineWidth', 1.8, 'Color', col2);

xlim([t(1) t(end)]);
ylim([0 1]);
xlabel('time');
ylabel('|z|');
title('Order parameter |z| of first/second half of oscillators');
legend([p1 p2], {'first half','second half'}, 'Location','best');
hold off;

%% ========= Helper: find most recent flag_complex_bundle.mat =========
function latest = findLatestBundle()
    latest = "";
    d = dir('flag_n200_p0.500*');
    bestTime = -inf;
    for k = 1:numel(d)
        if ~d(k).isdir, continue; end
        m = dir(fullfile(d(k).name, 'flag_complex_bundle.mat'));
        if ~isempty(m)
            [~, idx] = max([m.datenum]);
            if m(idx).datenum > bestTime
                bestTime = m(idx).datenum;
                latest = fullfile(d(k).name, m(idx).name);
            end
        end
    end
end
