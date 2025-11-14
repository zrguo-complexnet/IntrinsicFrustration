function [z1, z2] = simulate_oscillators(N, omega, t_max, dt, K)
    % simulate_oscillators
    % -------------------------------------------------------------
    % Simulates N phase oscillators with m-th order Kuramoto-type
    % coupling using a fourth-order Runge–Kutta scheme:
    %
    %     dθ_i/dt = ω_i − K R^n sin[n(θ_i − Θ)],
    %
    % where  Z = R e^{iΘ} = (1/N) Σ_j e^{iθ_j}.
    %
    % Inputs:
    %   N      – number of oscillators
    %   omega  – intrinsic frequencies (N×1 vector)
    %   t_max  – total simulation time
    %   dt     – time step
    %   K      – coupling strength
    %
    % Outputs:
    %   z1, z2 – complex order parameters of the two halves (N/2 each)
    % -------------------------------------------------------------

    n = 3;                         % Nonlinearity exponent
    t_steps = t_max / dt;          % Number of time steps

    % Initial phases: 10% at 2π/3, others at 0
    theta = (rand(N,1) < 0.1) * (2*pi/3);

    % Phase dynamics
    function dtheta = theta_dot(theta, omega, K, N, n)
        z     = mean(exp(1i*theta));        % Global order parameter
        R     = abs(z);
        Theta = angle(z);
        dtheta = omega - K*R^n .* sin(n*(theta - Theta));
    end

    % Time integration (RK4)
    for t = 1:t_steps
        k1 = dt * theta_dot(theta,              omega, K, N, n);
        k2 = dt * theta_dot(theta + 0.5*k1,     omega, K, N, n);
        k3 = dt * theta_dot(theta + 0.5*k2,     omega, K, N, n);
        k4 = dt * theta_dot(theta + k3,         omega, K, N, n);
        theta = theta + (k1 + 2*k2 + 2*k3 + k4)/6;
    end

    % Order parameters of the two subpopulations
    z1 = (2/N) * sum(exp(1i * theta(1:N/2)));
    z2 = (2/N) * sum(exp(1i * theta(N/2+1:N)));
end
