function [z1, z2] = Bsimulate_oscillators_K2(N, omega, t_max, dt, K, K2)
% BSIMULATE_OSCILLATORS_K2
% -------------------------------------------------------------------------
% Simulates N phase oscillators with mixed third- and second-harmonic
% mean-field Kuramoto-type coupling:
%
%   dθ_i/dt = ω_i
%             − K  R^3 sin[3(θ_i − Θ)]
%             − K1 R^2 sin[2(θ_i − Θ)],
%
% where Z = R e^{iΘ} = (1/N) Σ_j e^{iθ_j}.
%
% Inputs:
%   N      – number of oscillators
%   omega  – intrinsic frequencies (N×1 vector)
%   t_max  – total simulation time
%   dt     – time step
%   K      – third-order coupling strength
%   K2     – second-harmonic coupling strength
%
% Outputs:
%   z1, z2 – complex order parameters of the first and second halves
%            (each half contains N/2 oscillators)
% -------------------------------------------------------------------------

    n = 3;                         % Nonlinearity exponent for the highest-order term
    t_steps = t_max / dt;          % Number of time steps

    % Initial phases: symmetric two-group pattern
    unite1 = [0, 0, 0, 0, 2*pi/n, 0, 0, 0, 0, 0];
    theta1 = repmat(unite1, 1, N/20).';
    theta2 = flip(theta1);
    theta  = [theta1; theta2];

    % Phase dynamics with mixed 3rd and 2nd harmonics
    function dtheta = theta_dot(theta, omega, K, N, n)
        z     = (1/N) * sum(exp(1i * theta));   % Global order parameter
        R     = abs(z);
        phase = angle(z);
        dtheta = omega ...
               - (K  * R^n) .* sin(n * (theta - phase)) ...
               - (K2 * R^2) .* sin(2 * (theta - phase));
    end

    % Time integration via 4th-order Runge–Kutta
    for t = 1:t_steps
        k1 = dt * theta_dot(theta,            omega, K, N, n);
        k2 = dt * theta_dot(theta + 0.5*k1,   omega, K, N, n);
        k3 = dt * theta_dot(theta + 0.5*k2,   omega, K, N, n);
        k4 = dt * theta_dot(theta + k3,       omega, K, N, n);
        theta = theta + (k1 + 2*k2 + 2*k3 + k4) / 6;
    end

    % Order parameters for the two halves
    z1 = (2/N) * sum(exp(1i * theta(1:N/2)));       % First half
    z2 = (2/N) * sum(exp(1i * theta(N/2+1:N)));     % Second half
end
