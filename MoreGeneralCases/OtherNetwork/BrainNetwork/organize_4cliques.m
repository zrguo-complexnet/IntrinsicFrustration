function clique_matrix = organize_4cliques(four_cliques)
    n = size(four_cliques, 1); % 4-cliques的数量
    if n == 0
        error('没有找到4-cliques');
    end

    max_node = max(four_cliques(:)); % 查找最大节点编号
    clique_matrix = cell(max_node, 1); % 初始化单元格数组

    % 遍历每个4-clique，填充clique_matrix
    for idx = 1:n
        clique = four_cliques(idx, :);
        for node = clique
            clique_matrix{node} = [clique_matrix{node}, idx]; % 追加编号
        end
    end

    % 将单元格数组转换为矩阵
    max_cliques = cellfun(@numel, clique_matrix); % 每个节点的clique数量
    
    clique_matrix = cell2mat(cellfun(@(x) [x, NaN(1, max(max_cliques) - numel(x))], clique_matrix, 'UniformOutput', false)); 

    % 输出结果
    disp('Node | Connected 4-cliques');
    %disp(clique_matrix);
end