function [A, simplices, f] = randomFlagComplex(n, p, maxDim, rngSeed)
% RANDOMFLAGCOMPLEX
% ------------------------------------------------------------
% Generates an Erdős–Rényi graph G(n,p) and builds its flag
% (clique) complex up to dimension maxDim.
%
% INPUT
%   n        – number of vertices
%   p        – edge probability
%   maxDim   – maximal simplex dimension
%   rngSeed  – optional seed
%
% OUTPUT
%   A          – n×n adjacency matrix (0/1, symmetric)
%   simplices  – simplices{k} is the list of k-vertex cliques
%   f          – f(d+1) = number of d-dimensional simplices
% ------------------------------------------------------------

    if nargin < 3, maxDim = 3; end
    if nargin > 3 && ~isempty(rngSeed), rng(rngSeed); end

    % --- Generate G(n,p) ---
    A = triu(rand(n) < p, 1);
    A = logical(A + A.');
    A(1:n+1:end) = 0;

    % --- Clique expansion ---
    simplices = cell(maxDim+1,1);
    simplices{1} = (1:n).';
    f = zeros(maxDim+1,1);
    f(1) = n;

    for k = 2:maxDim+1
        prev = simplices{k-1};
        if isempty(prev), break; end

        list_k = [];
        for r = 1:size(prev,1)
            clique = prev(r,:);
            mask = true(n,1);
            for v = clique
                mask = mask & A(:,v);
            end
            cand = find(mask & ((1:n)' > clique(end)));
            if ~isempty(cand)
                list_k = [list_k; repmat(clique,numel(cand),1), cand(:)]; %#ok<AGROW>
            end
        end

        simplices{k} = list_k;
        f(k) = size(list_k,1);
        if isempty(list_k)
            simplices = simplices(1:k-1);
            f = f(1:k-1);
            break;
        end
    end
end
