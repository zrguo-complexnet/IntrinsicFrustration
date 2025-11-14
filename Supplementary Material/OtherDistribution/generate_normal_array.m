function result = generate_normal_array(N, center, sigma)
    % generate_normal_array
    % -------------------------------------------
    % Returns N samples from a normal distribution with
    % mean = center and standard deviation = sigma,
    % obtained by evaluating the inverse CDF (norminv)
    % at uniformly spaced quantiles.
    % -------------------------------------------

    result = norminv((1:(N+1)) / (N+1), center, sigma);
    result(N+1) = [];
end
