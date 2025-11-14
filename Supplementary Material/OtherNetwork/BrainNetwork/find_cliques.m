function [three_cliques, four_cliques] = find_cliques(A)
% FIND_CLIQUES  Brute-force enumeration of 3- and 4-cliques
% ----------------------------------------------------------
% Input:
%   A : n×n adjacency matrix (binary, symmetric, zero diagonal)
%
% Output:
%   three_cliques : all 3-cliques, each row = [i j k]
%   four_cliques  : all 4-cliques, each row = [i j k l]
%
% Note: O(n^3) and O(n^4) enumeration; intended for small graphs.
% ----------------------------------------------------------

    n = size(A,1);
    three_cliques = [];
    four_cliques  = [];

    % ----------- 3-cliques -----------
    for i = 1:n
        for j = i+1:n
            for k = j+1:n
                if A(i,j) && A(i,k) && A(j,k)
                    three_cliques = [three_cliques; i, j, k];
                end
            end
        end
    end

    % ----------- 4-cliques -----------
    for i = 1:n
        for j = i+1:n
            for k = j+1:n
                for l = k+1:n
                    if A(i,j) && A(i,k) && A(i,l) && ...
                       A(j,k) && A(j,l) && A(k,l)
                        four_cliques = [four_cliques; i, j, k, l];
                    end
                end
            end
        end
    end

    % Display (optional)
    disp('3-cliques:');
    disp(three_cliques);
    disp('4-cliques:');
    disp(four_cliques);
end
