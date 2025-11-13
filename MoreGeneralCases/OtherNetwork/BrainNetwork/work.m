% scale_free_network(6, 100, 3);
% random_network(100, 0.05);
% [three_cliques, four_cliques] = find_cliques(CIJ_resampled_average);
% clique_matrix = organize_4cliques(four_cliques); % 整理cliques

K2 = 15; % 耦合强度
K1 = 3;
num_steps = 5000; % 模拟步数
dt = 0.01; % 时间步长
kuramoto_simulation(clique_matrix, four_cliques, max_cliques, K1,K2, num_steps, dt,CIJ_resampled_average);
