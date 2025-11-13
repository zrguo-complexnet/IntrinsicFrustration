function result = generate_sorted_uniform_array(N, omega0, Delta)
    % 生成稳定排序的均匀分布数组
    % N: 数组的大小
    % omega0: 区间中心
    % Delta: 半宽度（结果落在 [omega0-Delta, omega0+Delta]）

    % 等距样本（已排序）
    U = linspace(0, 1, N);

    % 线性映射到目标区间
    result = omega0 + Delta * (2*U - 1);
end
