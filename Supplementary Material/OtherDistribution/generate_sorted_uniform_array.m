function result = generate_sorted_uniform_array(N, omega0, Delta)
    % generate_sorted_uniform_array
    % --------------------------------------------------------
    % Generates N sorted samples from a uniform distribution
    % supported on [omega0 - Delta, omega0 + Delta], obtained
    % by linearly mapping evenly spaced points in [0, 1].
    % --------------------------------------------------------

    U = linspace(0, 1, N);               % Sorted uniform samples in [0,1]
    result = omega0 + Delta*(2*U - 1);   % Map to target interval
end
