function [z1, z2] = Asimulate_oscillators_K1(N, omega, t_max, dt, K, K1)
% ASIMULATE_OSCILLATORS_K1
% -------------------------------------------------------------------------
% Simulates N phase oscillators with mixed first- and third-order
% mean-field Kuramoto-type coupling:
%
%   dθ_i/dt = ω_i
%             − K  R^3 sin[3(θ_i − Θ)]
%             − K1 R   sin[  (θ_i − Θ)],
%
% where Z = R e^{iΘ} = (1/N) Σ_j e^{iθ_j}.
%
% Inputs:
%   N      – number of oscillators
%   omega  – intrinsic frequencies (N×1 vector)
%   t_max  – total simulation time
%   dt     – time step
%   K      – third-order coupling strength
%   K1     – first-order coupling strength
%
% Outputs:
%   z1, z2 – complex order parameters of the first and second halves
%            (each half contains N/2 oscillators)
% -------------------------------------------------------------------------

    n = 3;                          % Nonlinearity exponent for the high-order term
    t_steps = t_max / dt;           % Number of time steps

    % Initial phases: symmetric two-group pattern
    unite1 = [0, 0, 0, 0, 2*pi/n, 0, 0, 0, 0, 0];
    theta1 = repmat(unite1, 1, N/20).';   % First half
    theta2 = flip(theta1);                % Second half (mirrored)
    theta  = [theta1; theta2];            % Full N×1 phase vector

    % Phase dynamics (mean-field)
    function dtheta = theta_dot(theta, omega, K, N, n)
        z     = (1/N) * sum(exp(1i * theta));   % Global order parameter
        R     = abs(z);
        phase = angle(z);
        dtheta = omega ...
               - (K  * R^n) .* sin(n * (theta - phase)) ...
               - (K1 * R   ) .* sin(    (theta - phase));
    end

    % Time integration via 4th-order Runge–Kutta
    for t = 1:t_steps
        k1 = dt * theta_dot(theta,              omega, K, N, n);
        k2 = dt * theta_dot(theta + 0.5 * k1,   omega, K, N, n);
        k3 = dt * theta_dot(theta + 0.5 * k2,   omega, K, N, n);
        k4 = dt * theta_dot(theta + k3,         omega, K, N, n);
        theta = theta + (k1 + 2*k2 + 2*k3 + k4) / 6;
    end

    % Order parameters for the two halves
    z1 = (2/N) * sum(exp(1i * theta(1:N/2)));         % First half
    z2 = (2/N) * sum(exp(1i * theta(N/2+1:N)));       % Second half
end
