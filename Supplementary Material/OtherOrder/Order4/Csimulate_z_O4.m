clear;
clc;

N      = 1000;    % Number of oscillators
t_max  = 20;      % Total simulation time
dt     = 1e-4;    % Time step

if isempty(gcp('nocreate'))
    parpool(30);  % Start a parallel pool if none exists
end

% Parameter ranges
K_values      = linspace(0, 900, 100);    % K from 0 to 900 (100 values)
omega0_values = linspace(0, 1200, 100);   % Frequency offset from 0 to 1200

z_diff_matrix = zeros(length(K_values), length(omega0_values));  % (optional) |z1 - z2| storage

% Base Lorentzian-shaped sample, centered at 0
lorentz = generate_sorted_lorentz_array(N/2, 0, 1)';

SL = length(omega0_values);

% Parameter sweep over (K, omega0)
parfor i = 1:length(K_values)
    for j = 1:SL
        K = K_values(i);  % Current 4th-order coupling strength

        % Local copy of the omega0 grid (kept to minimally modify code)
        omega0_values_local = linspace(0, 1200, 100); %#ok<NASGU>
        omega0 = omega0_values(j);  % Current frequency separation

        % Construct intrinsic frequency vector (bimodal, symmetric)
        omega = [ lorentz + omega0; ...
                 -lorentz - omega0 ];

        % Simulate oscillator system with pure 4th-order coupling
        [z1, z2] = simulate_oscillators1(N, omega, t_max, dt, K);

        % Store complex order parameters of the two halves
        % z_diff_matrix(i, j) = abs(abs(z2) - abs(z1));  % optional magnitude difference
        z_1(i, j) = z1;
        z_2(i, j) = z2;
    end
end

mm1 = z_1;
mm2 = z_2;

% Save results
save('z1Order4Sim910D1_9001200.mat', 'mm1');
save('z2Order4Sim910D1_9001200.mat', 'mm2');

delete(gcp('nocreate'));   % Shut down the parallel pool
