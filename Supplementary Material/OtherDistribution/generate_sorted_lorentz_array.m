function result = generate_sorted_lorentz_array(N, omega0, Delta)
    % generate_sorted_lorentz_array
    % --------------------------------------------------------
    % Generates N deterministically sorted samples from a
    % Lorentzian (Cauchy) distribution with location omega0
    % and scale Delta, using the inverse CDF method evaluated
    % at uniformly spaced quantiles.
    % --------------------------------------------------------

    U = linspace(0.01, 0.99, N);                 % Uniform quantiles
    result = omega0 + Delta * tan(pi*(U - 0.5)); % Inverse CDF of Lorentzian
end
