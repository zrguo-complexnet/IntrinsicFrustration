function [three_cliques, four_cliques] = find_cliques(A)
    % A 是邻接矩阵
    n = size(A, 1); % 节点数量
    three_cliques = []; % 存储3-clique
    four_cliques = [];  % 存储4-clique

    % 查找3-clique
    for i = 1:n
        for j = i+1:n
            for k = j+1:n
                if A(i,j) && A(i,k) && A(j,k) % 检查是否是完全连接
                    three_cliques = [three_cliques; i, j, k];
                end
            end
        end
    end

    % 查找4-clique
    for i = 1:n
        for j = i+1:n
            for k = j+1:n
                for l = k+1:n
                    if A(i,j) && A(i,k) && A(i,l) && A(j,k) && A(j,l) && A(k,l) % 检查是否是完全连接
                        four_cliques = [four_cliques; i, j, k, l];
                    end
                end
            end
        end
    end

    % 输出结果
    disp('3-cliques:');
    disp(three_cliques);
    disp('4-cliques:');
    disp(four_cliques);
end
