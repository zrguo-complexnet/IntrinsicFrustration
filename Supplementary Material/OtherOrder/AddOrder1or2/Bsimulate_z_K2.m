clear;
clc;

N      = 1000;   % Number of oscillators
t_max  = 50;     % Total simulation time
dt     = 0.001;  % Time step
omega0 = 50;     % Frequency shift between the two groups

if isempty(gcp('nocreate'))
    parpool(30); % Start a parallel pool if none exists
end

% Parameter ranges
K_values  = linspace(0, 160, 100);  % K from 0 to 160 (100 values)
K1_values = linspace(0, 10, 100);   % K1 from 0 to 10 (100 values)

z_diff_matrix = zeros(length(K_values), length(K1_values));  % (optional) |z1 - z2| storage

% Generate a symmetric Lorentzian-like distribution for intrinsic frequencies
lorentz = generate_sorted_lorentz_array(N/2, 0, 1)';

SL = length(K1_values);

% Construct intrinsic frequencies: +omega0 for the first half, -omega0 for the second
omega = [ lorentz + omega0; ...
         -lorentz - omega0 ];

% Sweep over (K, K1)
parfor i = 1:length(K_values)
    for j = 1:SL
        K = K_values(i);   % Current third-order coupling strength

        % Local copy of K1 grid inside the loop (kept for compatibility;
        % note: here it is redefined from 0 to 5)
        K1_values_local = linspace(0, 5, 100);  %#ok<NASGU>
        K1 = K1_values(j);  % Current first-order coupling strength

        % Simulate the oscillator system
        [z1, z2] = Asimulate_oscillators_K1(N, omega, t_max, dt, K, K1);

        % Store complex order parameters
        % z_diff_matrix(i, j) = abs(abs(z2) - abs(z1));  % optional measure
        z_1(i, j) = z1;
        z_2(i, j) = z2;
    end
end

mm1 = z_1;
mm2 = z_2;

% Save data to file
save('z1Sim910D1_K2_16010.mat', 'mm1');
save('z2Sim910D1_K2_16010.mat', 'mm2');

delete(gcp('nocreate'));  % Shut down the parallel pool
