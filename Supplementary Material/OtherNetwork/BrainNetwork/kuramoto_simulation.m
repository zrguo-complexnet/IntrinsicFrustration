function kuramoto_simulation(clique_matrix, four_cliques, max_cliques, ...
                             K1, K2, num_steps, dt, CIJ_resampled_average)
% KURAMOTO_SIMULATION
% -------------------------------------------------------------------------
% Kuramoto dynamics with:
%   - Pairwise coupling on CIJ_resampled_average (strength K1)
%   - 4-body coupling on hyperedges given by four_cliques (strength K2),
%     normalized by node-wise 4-clique counts max_cliques.
%
% Inputs:
%   clique_matrix           – N×M index matrix linking nodes to 4-cliques
%   four_cliques            – (#4-cliques)×4 array of node indices
%   max_cliques             – N×1, number of 4-cliques incident to each node
%   K1, K2                  – pairwise and 4-body coupling strengths
%   num_steps               – number of time steps
%   dt                      – time step
%   CIJ_resampled_average   – N×N binary adjacency for pairwise coupling
% -------------------------------------------------------------------------

    R1_values = zeros(num_steps, 1);
    R2_values = zeros(num_steps, 1);

    N = size(clique_matrix, 1);   % number of oscillators

    % Initial phases: repeated zero pattern
    unite1 = zeros(1,20);
    theta1 = repmat(unite1, 1, (N+2)/40);
    theta1((N+2)/40) = [];
    theta = [theta1'; theta1'];   % two identical halves

    % Bimodal intrinsic frequencies
    omega = generate_bimodal_distribution(N);

    % ---------------------- Time evolution ----------------------
    for step = 1:num_steps
        dtheta = zeros(N, 1);

        % ----- 4-body contributions via four_cliques -----
        for i = 1:N
            d      = 0;
            dtheta1 = 0;

            if max_cliques(i) > 0
                for j = clique_matrix(i,:)
                    if isnan(j)
                        continue;
                    end
                    cliques = four_cliques(j, :);
                    dtheta(i) = dtheta(i) + ...
                        sin(sum(theta(cliques)) - 4*theta(i));
                end
                dtheta(i) = omega(i) + (K2 / max_cliques(i)) * dtheta(i);
            end

            % ----- Pairwise Kuramoto term on CIJ_resampled_average -----
            for j = 1:N
                if CIJ_resampled_average(i,j)
                    d      = d + 1;
                    dtheta1 = dtheta1 + sin(theta(j) - theta(i));
                end
            end
            if d
                dtheta(i) = dtheta(i) + (K1/d)*dtheta1;
            end
        end

        % Order parameters for the two halves
        R  = abs((1/N) * sum(exp(1i * theta)));          %#ok<NASGU>
        R1 = abs((2/N) * sum(exp(1i * theta(1:N/2))));
        R2 = abs((2/N) * sum(exp(1i * theta(N/2+1:N))));

        R1_values(step) = R1;
        R2_values(step) = R2;

        % Phase update
        theta = theta + dtheta * dt;
        theta = mod(theta, 2*pi);
    end

    % ---------------------- Plot R1, R2 vs time ----------------------
    figure;
    plot(1:num_steps, R1_values, 'b-', 'LineWidth', 1.5);
    hold on;
    plot(1:num_steps, R2_values, 'r-', 'LineWidth', 1.5);
    xlabel('Time step');
    ylabel('R');
    title('R_1 and R_2 over time');
    legend('R_1', 'R_2');
    grid on;
end

function omega = generate_bimodal_distribution(N)
% GENERATE_BIMODAL_DISTRIBUTION
% Generates a bimodal normal distribution: N/2 around -5 and N/2 around +5.

    mu1 = -5; sigma1 = 0.1;
    mu2 =  5; sigma2 = 0.1;
    omega = [normrnd(mu1, sigma1, [N/2, 1]); ...
             normrnd(mu2, sigma2, [N/2, 1])];
    % omega = omega(randperm(N));  % optional shuffling
end
