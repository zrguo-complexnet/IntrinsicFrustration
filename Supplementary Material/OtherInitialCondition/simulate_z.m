clear; clc;

N      = 1000;          % Number of oscillators
t_max  = 20;            % Total simulation time
dt     = 1e-4;          % Time step

% Start a parallel pool if none exists
if isempty(gcp('nocreate'))
    parpool(30);
end

% ---------------------------- Parameter grids -----------------------------
K_values      = linspace(0, 200, 100);    % Coupling strength K
omega0_values = linspace(0, 150, 100);    % Frequency shift parameter
SL = length(omega0_values);

% Base distribution for natural frequencies (choose distribution as needed)
lorentz = generate_normal_array(N/2, 0, 1)';   % Gaussian samples (example)

% Preallocate storage for the two order parameters
z_1 = zeros(length(K_values), SL);
z_2 = zeros(length(K_values), SL);

% --------------------------- Two-parameter sweep --------------------------
parfor i = 1:length(K_values)
    K = K_values(i);
    for j = 1:SL

        omega0 = omega0_values(j);

        % Construct natural frequencies:
        %   first half  :  lorentz + omega0
        %   second half : -lorentz - omega0
        omega = [lorentz + omega0;
                -lorentz - omega0];

        % Simulate dynamics and compute order parameters (complex)
        [z1, z2] = simulate_oscillators(N, omega, t_max, dt, K);

        z_1(i, j) = z1;
        z_2(i, j) = z2;
    end
end

mm1 = z_1;
mm2 = z_2;

% ------------------------------- Save data --------------------------------
save('z1Sim910D1_150200_.mat', 'mm1');
save('z2Sim910D1_150200_.mat', 'mm2');

% Close parallel pool
delete(gcp('nocreate'));
