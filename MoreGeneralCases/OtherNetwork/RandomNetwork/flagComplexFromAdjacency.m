function [A, simplices, f] = flagComplexFromAdjacency(A, maxDim)
% FLAGCOMPLEXFROMADJACENCY  从给定邻接矩阵 A（n×n、0/1、无向无自环）
%                           构造其随机团（flag/clique）复形到最高维 maxDim。
% 输出格式与 randomFlagComplex 一致：
%   simplices{1} = (1:n)'，simplices{2}=所有边，simplices{3}=所有三角形，...

    if nargin<2 || isempty(maxDim), maxDim = 3; end
    n = size(A,1);
    A = logical(A);
    A = triu(A,1) | triu(A,1).';  % 确保对称无向
    A(1:n+1:end) = 0;             % 无自环

    simplices = cell(maxDim+1,1);
    simplices{1} = (1:n).';
    f = zeros(maxDim+1,1);  f(1) = n;

    for k = 2:(maxDim+1)
        prev = simplices{k-1};
        if isempty(prev), break; end

        list_k = [];
        for r = 1:size(prev,1)
            clique = prev(r,:);         % 已知 (k-1)-团
            mask = true(n,1);
            for t = 1:numel(clique)
                mask = mask & A(:,clique(t));
            end
            cand = find(mask & ((1:n)' > clique(end)));
            if ~isempty(cand)
                new_rows = [repmat(clique, numel(cand), 1), cand(:)];
                list_k = [list_k; new_rows]; %#ok<AGROW>
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
