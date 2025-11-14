function [clique_matrix, max_cliques] = organize_4cliques(four_cliques)
    n = size(four_cliques, 1); % number of 4-cliques
    if n == 0
        error('No 4-cliques found.');
    end

    max_node = max(four_cliques(:));         % largest node index
    clique_cell = cell(max_node, 1);         % temporary cell array

    % For each 4-clique, append its index to the incident nodes
    for idx = 1:n
        clique = four_cliques(idx, :);
        for node = clique
            clique_cell{node} = [clique_cell{node}, idx];
        end
    end

    % Number of 4-cliques incident to each node
    max_cliques = cellfun(@numel, clique_cell);

    % Convert cell array to an N×max_deg matrix, NaN-padded
    max_deg = max(max_cliques);
    clique_matrix = cell2mat(cellfun( ...
        @(x) [x, NaN(1, max_deg - numel(x))], ...
        clique_cell, 'UniformOutput', false));

end
