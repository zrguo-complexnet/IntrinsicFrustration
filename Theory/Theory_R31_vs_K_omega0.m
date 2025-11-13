% Define parameters
gamma1 = 0.9;
gamma2 = 0.1;
gamma3 = 0;     % Proportions of oscillators in the three synchronized regions
dt = 0.0001;    % Time step, should be as small as possible
N = 100000;     % Total number of time steps
Delta = 1;

load('R1_910.mat');
load('R3_910.mat');     % Data from the file "R1_R3.m"
R1_data = solutions;    % R1 data points
R3_data = solutions3;   % R3 data points

% Construct a functional relationship between R_{1,p} and R_{3,p} within each synchronized region
R1_R3_interp = @(R3) interp1(R3_data, R1_data, R3, 'spline');

% Enable parallel pool if not already active
if isempty(gcp('nocreate'))
    parpool(24);
end

% Define the ranges of K and omega0
K_values = linspace(0, 200, 100);
omega0_values = linspace(0, 150, 100);

% Initialize result matrix
results = zeros(length(K_values), length(omega0_values));
results2 = zeros(length(K_values), length(omega0_values));

SO = length(omega0_values);

% Loop over K and omega0 values
parfor i = 1:length(K_values)
    for j = 1:SO
        K = K_values(i);
        omega0 = omega0_values(j);
        
        % Initial conditions
        R31 = 1;
        Theta31 = 0;
        R32 = 1;
        Theta32 = 0;
        
        % Compute the partial order parameter z_1 from three synchronized regions,
        % assuming Theta_{1p} ≈ Theta_{3p} / 3
        z11 = R1_R3_interp(R31) * (gamma1 * exp(1i * Theta31 / 3) + gamma2 * exp(1i * (Theta31 / 3 + 2*pi/3)) + gamma3 * exp(1i * (Theta31 / 3 + 4*pi/3)));
        z12 = R1_R3_interp(R32) * (gamma1 * exp(1i * Theta32 / 3) + gamma2 * exp(1i * (Theta32 / 3 + 2*pi/3)) + gamma3 * exp(1i * (Theta32 / 3 + 4*pi/3)));
        z1 = 0.5 * (z11 + z12);  % Total order parameter is the average of the two subpopulations
        R1 = abs(z1);
        Theta1 = angle(z1);
        
        % Simulation loop
        for k = 1:N
            % Compute derivatives
            dR31 = (3*K/2) * R1^3 * cos(3*Theta1 - Theta31) * (1 - R31^2) - 3*R31*Delta;
            dTheta31 = (3*K/2) * R1^3 * sin(3*Theta1 - Theta31) * (R31 + 1/R31) - 3*omega0;
            dR32 = (3*K/2) * R1^3 * cos(3*Theta1 - Theta32) * (1 - R32^2) - 3*R32*Delta;
            dTheta32 = (3*K/2) * R1^3 * sin(3*Theta1 - Theta32) * (R32 + 1/R32) + 3*omega0;
            
            % Update variables
            R31 = R31 + dR31 * dt;
            Theta31 = Theta31 + dTheta31 * dt;
            R32 = R32 + dR32 * dt;
            Theta32 = Theta32 + dTheta32 * dt;

            % Recalculate z1, R1, and Theta1
            z11 = R1_R3_interp(R31) * (gamma1 * exp(1i * Theta31 / 3) + gamma2 * exp(1i * (Theta31 / 3 + 2*pi/3)) + gamma3 * exp(1i * (Theta31 / 3 + 4*pi/3)));
            z12 = R1_R3_interp(R32) * (gamma1 * exp(1i * Theta32 / 3) + gamma2 * exp(1i * (Theta32 / 3 + 2*pi/3)) + gamma3 * exp(1i * (Theta32 / 3 + 4*pi/3)));
            z1 = 0.5 * (z11 + z12);
            R1 = abs(z1);
            Theta1 = angle(z1);
        end
        
        % Store result
        results(i, j) = abs(R31);
        results2(i, j) = abs(R32);
    end
end

save('z1_D1_TR910_200_150.mat', 'results');  % Save data to file
save('z2_D1_TR910_200_150.mat', 'results');  % Save data to file
delete(gcp('nocreate'))

