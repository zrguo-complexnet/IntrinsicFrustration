function kuramoto_simulation(clique_matrix, four_cliques, max_cliques, K1,K2, num_steps, dt,CIJ_resampled_average)
    % 设置参数
    R1_values = zeros(num_steps, 1);
    R2_values = zeros(num_steps, 1);
    N = size(clique_matrix, 1); % 节点数量
     unite1 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
% %     unite2 = [0,0,0,0,0,2*pi/3,0,0,0,0];
     theta1 = repmat(unite1, 1, (N+2)/40);
     theta1((N+2)/40) = [];
%      theta2 = rand(N/2,1)*2*pi;
% %     theta2 = repelem(unite2, 1, N/20)';
     theta = [theta1'; theta1'] ;
%    theta = zeros(N, 1); % 初始相位
    omega = generate_bimodal_distribution(N); % 双峰正态分布

    % 动力学模拟
%     figure;
    for step = 1:num_steps
        dtheta = zeros(N, 1); % 每个节点的相位变化率
        
        for i = 1:N
            d=0;
            dtheta1=0;
            if max_cliques(i) > 0
                for j = clique_matrix(i,:)
                    if isnan(j)
                        continue; % 跳过未定义的点
                    end
                    % 查找与节点i连接的所有4-clique
                    cliques = four_cliques(j, :);
                    % 对每个连接的4-clique进行贡献计算
                    dtheta(i) = dtheta(i) + sin(sum(theta(cliques)) - 4 * theta(i));
                end
                dtheta(i) = omega(i)+(K2 / max_cliques(i)) * dtheta(i);
            end
            for j=1:N
                if CIJ_resampled_average(i,j)
                    d=d+1;
                    dtheta1=dtheta1+sin(theta(j)-theta(i));
                end
            end
            if d 
                dtheta(i)=dtheta(i)+(K1/d)*dtheta1;
            end
        end
        R = abs((1/N) * sum(exp(1i * theta)));
        R1 = abs((2/N) * sum(exp(1i * theta(1:N/2))));
        R2 = abs((2/N) * sum(exp(1i * theta(N/2+1:N))));
        
        R1_values(step) = R1;
        R2_values(step) = R2;
        
        theta = theta + dtheta * dt; % 更新相位
        theta = mod(theta,2*pi);
        % 绘制振子的相位在圆上
%         clf; % 清除当前图形
%         hold on;
%         % 计算振子在圆上的坐标
%         x = omega;
%         y = theta;
%         plot(x, y, 'o', 'MarkerSize', 8); % 绘制振子
%         axis equal; % 设置坐标轴等比例
%         xlim([-6, 6]);
%         ylim([0, 2*pi]);
%         title(sprintf('Kuramoto Model Dynamics - Step %d - R %f\n R1 %f - R2 %f', step, R, R1, R2));
%         xlabel('X');
%         ylabel('Y');
%         grid on;
% %         pause(0.1); % 控制动画速度
    end
    
    % 绘制R1和R2的随时间变化图
    figure; % 创建新的图形窗口
    plot(1:num_steps, R1_values, 'b-', 'LineWidth', 1.5);
    hold on;
    plot(1:num_steps, R2_values, 'r-', 'LineWidth', 1.5);
    xlabel('Time Step');
    ylabel('R Value');
    title('R1 and R2 Dynamics Over Time');
    legend('R1', 'R2');
    grid on;
end

function omega = generate_bimodal_distribution(N)
    % 生成双峰正态分布
    mu1 = -5; sigma1 = 0.1; % 第一个峰的均值和标准差
    mu2 = 5; sigma2 = 0.1; % 第二个峰的均值和标准差
    omega = [normrnd(mu1, sigma1, [N/2, 1]); normrnd(mu2, sigma2, [N/2, 1])];
   % omega = omega(randperm(N)); % 打乱顺序
end
