function [A, simplices, f] = randomFlagComplex(n, p, maxDim, rngSeed)
% RANDOMFLAGCOMPLEX  生成 G(n,p) 并构造其随机团（flag/clique）复形
% 输入:
%   n        - 节点数
%   p        - 边独立出现概率
%   maxDim   - 复形的最高维 (>=1 表示至少到边；2含三角形；3含四面体; …)
%   rngSeed  - (可选) 随机种子，便于复现实验
% 输出:
%   A          - n×n 的邻接矩阵（0/1，对称、无自环）
%   simplices  - 单元格数组，simplices{k} 是所有 (k)-节点团（即 (k-1) 维单纯形）
%                例如 simplices{1} = 所有顶点 (n×1)；
%                     simplices{2} = 所有边 (m2×2)；
%                     simplices{3} = 所有三角形 (m3×3)；依此类推
%   f          - f-向量，f(d+1) = d维单纯形个数（顶点记作0维）
%
% 复杂度说明:
%   团枚举本质是指数级；请用合适的 maxDim（常见到 2 或 3）。
%
% 用法示例（脚本里）:
%   [A,S,f] = randomFlagComplex(200, 0.05, 3, 42);
%   fprintf('顶点=%d, 边=%d, 三角形=%d, 四面体=%d\n', f(1), f(2), f(3), f(4));

    if nargin < 3 || isempty(maxDim), maxDim = 3; end
    if nargin >= 4 && ~isempty(rngSeed), rng(rngSeed); end

    %--- 1) 生成 G(n,p) ---
    A = triu(rand(n) < p, 1);
    A = logical(A + A.');        % 对称化；逻辑矩阵更快
    A(1:n+1:end) = 0;            % 去自环

    %--- 2) 自底向上枚举团（直到 maxDim） ---
    simplices = cell(maxDim+1,1);
    simplices{1} = (1:n).';      % 0维：所有顶点（作为 1 列）
    f = zeros(maxDim+1,1);
    f(1) = n;

    % 逐维扩展：把 (k-1)-团扩成 k-团（对应 (k-1) 维单纯形到 k-1 维的）
    for k = 2:(maxDim+1)
        prev = simplices{k-1};     % 每行一个已知团（大小 k-1）
        if isempty(prev), break; end

        list_k = [];               % 将收集所有 k-团（k 个顶点）
        list_k_capacity = 0;

        % 为了去重：仅向“比最后一个顶点编号更大的顶点”扩展
        for r = 1:size(prev,1)
            clique = prev(r,:);                 % 1×(k-1)
            % 计算该团所有顶点的公共邻居
            mask = true(n,1);
            for t = 1:numel(clique)
                mask = mask & A(:,clique(t));
            end
            % 只允许扩展到编号更大的顶点
            cand = find(mask & ((1:n)' > clique(end)));
            if ~isempty(cand)
                new_rows = [repmat(clique, numel(cand), 1), cand(:)];
                % 预分配/累加
                if isempty(list_k)
                    list_k = new_rows;
                    list_k_capacity = size(new_rows,1);
                else
                    list_k = [list_k; new_rows]; %#ok<AGROW>
                    list_k_capacity = list_k_capacity + size(new_rows,1);
                end
            end
        end

        % 保存本维度的所有单纯形（k 个顶点 ⇒ (k-1) 维单纯形）
        simplices{k} = list_k;
        f(k) = size(list_k,1);
        if isempty(list_k)
            % 没有更高维的团了，提前结束
            simplices = simplices(1:k-1);
            f = f(1:k-1);
            break;
        end
    end
end
