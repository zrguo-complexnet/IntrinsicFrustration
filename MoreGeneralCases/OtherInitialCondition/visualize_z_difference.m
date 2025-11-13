clear;
clc;

N = 1000;  % 振子的数量
t_max = 20;
dt = 0.0001;

if isempty(gcp('nocreate'))
    parpool(30);
end

% 初始化参数范围
K_values = linspace(0,200, 100);  % K从0到100，取100个值
omega0_values = linspace(0, 150, 100);  % omega0从0到100，取100个值
z_diff_matrix = zeros(length(K_values), length(omega0_values));  % 用来存储z1-z2差异
lorentz = generate_normal_array(N/2, 0, 1)'; %在这里可以选择不同的分布类型

SL=length(omega0_values);
% 遍历每一组K和标准差
parfor i = 1:length(K_values)
    for j = 1:SL
        K = K_values(i);  % 当前K值
            omega0_values = linspace(0, 150, 100);  % 标准差从0.1到2，取20个值
        omega0 = omega0_values(j);  % 当前距离

        % 生成自然频率数组omega，均值为0，标准差为std_dev
        omega = [lorentz+omega0;-lorentz-omega0] ;


        % 模拟振子系统
        [z1,z2] = simulate_oscillators(N, omega, t_max, dt, K);

        % 计算|z1 - z2|
%         z_diff_matrix(i, j) = abs(abs(z2)-abs(z1));
        z_1(i, j) = z1;
        z_2(i, j) = z2;
    end
end
mm1 =z_1;
mm2 =z_2;

save('z1Sim910D1_150200_.mat', 'mm1'); % 保存数据到文件 
save('z2Sim910D1_150200_.mat', 'mm2'); % 保存数据到文件 

delete(gcp('nocreate'))

