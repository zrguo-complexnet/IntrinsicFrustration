function [z1, z2] = simulate_oscillators(N, omega, t_max, dt, K)
    % simulate_oscillators
    % -------------------------------------------------------------
    % Simulates a population of N phase oscillators governed by an
    % n-th order sinusoidal coupling, using a fourth-order Runge–
    % Kutta time integrator. The population is divided into two
    % equal groups for computing the order parameters z1 and z2.
    %
    % Inputs:
    %   N      – total number of oscillators (even)
    %   omega  – intrinsic frequencies (vector of length N)
    %   t_max  – total simulation time
    %   dt     – time step
    %   K      – coupling strength
    %
    % Outputs:
    %   z1     – complex order parameter of the first half
    %   z2     – complex order parameter of the second half
    % -------------------------------------------------------------

    n = 3;                         % exponent in the nonlinear coupling
    t_steps = t_max / dt;          % number of time steps

    % Initial phases: structured pattern repeated across oscillators
    base = [0, 0, 0, 0, 2*pi/n, 0, 0, 0, 0, 0];
    theta1 = repmat(base, 1, N/20)';      % first group
    theta2 = flip(theta1);                % mirrored second group
    theta = [theta1; theta2];             % full initial condition

    % Phase dynamics
    function dtheta = theta_dot(theta, omega, K, N, n)
        z = (1/N) * sum(exp(1i*theta));       % global order parameter
        R = abs(z);
        phi = angle(z);
        dtheta = omega - K*R^n .* sin(n*(theta - phi));
    end

    % Time integration (RK4)
    for t = 1:t_steps
        k1 = dt * theta_dot(theta, omega, K, N, n);
        k2 = dt * theta_dot(theta + 0.5*k1, omega, K, N, n);
        k3 = dt * theta_dot(theta + 0.5*k2, omega, K, N, n);
        k4 = dt * theta_dot(theta + k3,    omega, K, N, n);
        theta = theta + (k1 + 2*k2 + 2*k3 + k4)/6;
    end

    % Order parameters of the two subpopulations
    z1 = (2/N) * sum(exp(1i*theta(1:N/2)));
    z2 = (2/N) * sum(exp(1i*theta(N/2+1:N)));
end
