function [z1,z2] = simulate_oscillators(N, omega, t_max, dt, K)
    % 参数设定
    n = 3;  % 非线性项指数
    t_steps = t_max / dt;  % 时间步数
%     unite1 = [0,0,0,0,2*pi/n,0,0,0,0,0,0,0,0,0,2*pi/n,0,0,0,0,0];
%     unite1 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
    % rd = rand(N/2,1)*2*pi/3;
    % unite1 = [0,0,0,0,2*pi/n,0,0,0,0,0];
    % theta1 = repmat(unite1, 1, N/20)';
%     theta1 =(rand(N/2, 1) < 0.1)*2/3*pi;
    % theta2=flip(theta1);
%     theta2 = zeros(N,1);
    % theta = [theta1; theta2] ;
%     theta = zeros(N, 1);  % 振子的初始相位为0
    theta =rand(N, 1)*0+(rand(N, 1) < 0.1)*2/3*pi;


    % 4阶Runge-Kutta方法定义
    function dtheta = theta_dot(theta, omega, K, N, n)
        z = (1/N) * sum(exp(1i * theta));  % 计算集体参数z
        R = abs(z);  % R是z的模长
        phase = angle(z);  % 集体相位
        dtheta = omega - (K * R^n) * sin(n * (theta - phase));  % 相位变化率
    end

    % 开始时间演化
    for t = 1:t_steps
        % Runge-Kutta 4阶更新过程
        k1 = dt * theta_dot(theta, omega, K, N, n);
        k2 = dt * theta_dot(theta + 0.5 * k1, omega, K, N, n);
        k3 = dt * theta_dot(theta + 0.5 * k2, omega, K, N, n);
        k4 = dt * theta_dot(theta + k3, omega, K, N, n);
        theta = theta + (k1 + 2*k2 + 2*k3 + k4) / 6;  % 更新相位
        %theta = mod(theta+pi/3,2*pi/3)-pi/3;
    end

    % 计算z1和z2
    z1 = (2/N) * sum(exp(1i * theta(1:N/2)));  % 前半部分振子
    z2 = (2/N) * sum(exp(1i * theta(N/2+1:N)));  % 前半部分振子
end
