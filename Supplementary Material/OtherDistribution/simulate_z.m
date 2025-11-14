clear; clc;

N     = 1000;          % Total number of oscillators
t_max = 20;            % Simulation horizon
dt    = 1e-4;          % Time step

% Open a parallel pool if none exists
if isempty(gcp('nocreate'))
    parpool;
end

% Parameter grids for the two-dimensional sweep
K_values     = linspace(0, 100, 100);   % Coupling strength K
omega0_values = linspace(0, 100, 100);  % Frequency shift parameter
z_diff_matrix = zeros(length(K_values), length(omega0_values));

% Base distribution for intrinsic frequencies (choose distribution here)
lorentz = generate_normal_array(N/2, 0, 1)';   % Example: Gaussian samples

SL = length(omega0_values);

% --------------------------------------------------------------
% Two-parameter sweep over (K, omega0)
% The intrinsic frequencies of the two subpopulations are
%    omega(1:N/2)     =  lorentz + omega0
%    omega(N/2+1:N)   = -lorentz - omega0
% --------------------------------------------------------------
parfor i = 1:length(K_values)
    for j = 1:SL

        K = K_values(i);                 % Current coupling strength
        omega0 = omega0_values(j);       % Current frequency shift

        % Construct the full intrinsic frequency vector
        omega = [lorentz + omega0; 
                -lorentz - omega0];

        % Run oscillator dynamics and compute order parameters
        [z1, z2] = simulate_oscillators(N, omega, t_max, dt, K);

        % Store full complex-valued order parameters
        z_1(i, j) = z1;
        z_2(i, j) = z2;
    end
end

mm1 = z_1;
mm2 = z_2;

% Save simulation results
save('z1GauSim910D1_100100.mat', 'mm1');
save('z2GauSim910D1_100100.mat', 'mm2');

% Close parallel pool
delete(gcp('nocreate'));
