clear;
clc;

N = 1000;  % 振子的数量
t_max = 50;
dt = 0.001;
omega0 = 50;

if isempty(gcp('nocreate'))
    parpool(30);
end

% 初始化参数范围
K_values = linspace(0,160, 100);  % K从0到5，取20个值
K1_values = linspace(0, 10, 100);  % 标准差从0.1到2，取20个值
z_diff_matrix = zeros(length(K_values), length(K1_values));  % 用来存储z1-z2差异
lorentz = generate_sorted_lorentz_array(N/2, 0, 1)';

SL=length(K1_values);
% 生成自然频率数组omega，均值为0，标准差为std_dev
omega = [lorentz+omega0;-lorentz-omega0] ;

% 遍历每一组K和标准差
parfor i = 1:length(K_values)
    for j = 1:SL
        K = K_values(i);  % 当前K值
        K1_values = linspace(0, 5, 100);  % 标准差从0.1到2，取20个值
        K1 = K1_values(j);  % 当前距离

        % 模拟振子系统
        [z1,z2] = Asimulate_oscillators_K1(N, omega, t_max, dt, K, K1);

        % 计算|z1 - z2|
%         z_diff_matrix(i, j) = abs(abs(z2)-abs(z1));
        z_1(i, j) = z1;
        z_2(i, j) = z2;
    end
end
mm1 =z_1;
mm2 =z_2;

save('z1Sim910D1_K2_16010.mat', 'mm1'); % 保存数据到文件 
save('z2Sim910D1_K2_16010.mat', 'mm2'); % 保存数据到文件 

delete(gcp('nocreate'))

