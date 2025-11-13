function result = generate_sorted_lorentz_array(N, omega0, Delta)

    % 生成稳定排序的Lorentzian分布数组
    % N: 数组的大小
    % omega0: Lorentzian分布的位置参数
    % Delta: Lorentzian分布的尺度参数
    
    % 生成从0到1的均匀分布的随机数
    U = linspace(0.01, 0.99, N);
    
    % 使用逆CDF公式计算Lorentzian分布的值
    result = omega0 + Delta * tan(pi * (U - 0.5));

end