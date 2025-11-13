% ==========================================================
% run_highorder_kuramoto_on_flag.m
% 在已构建好的随机团（flag）复形上，进行至 6 体的高阶 Kuramoto 演化
% 规则：
%   - 前半振子： +omega0 + sigma*rand
%   - 后半振子： -omega0 - sigma*rand
%   - r 阶相互作用：每个节点的贡献除以该节点的 r 阶度（连接的 r 阶超边数）
%   - 初相：全 0
% ==========================================================
clear; clc;

%% ========= 1) 选择要用的网络（.mat）与仿真参数 =========
bundleMat = "";  % 留空自动找最近一次输出的 flag_* 目录下的 flag_complex_bundle.mat
if bundleMat == "", bundleMat = string(findLatestBundle()); end
if bundleMat == "", error('未找到 flag_* 目录下的 flag_complex_bundle.mat'); end
Sdata = load(bundleMat);
fprintf('Loaded: %s\n', bundleMat);

% --- Kuramoto 参数 ---
omega0 = 5;      % 你的 +/− 基频
sigma  = 0.1;      % 频率扰动强度（与 rand 相乘）
Kvec   = zeros(6,1);
%      Kvec(2..6) 分别对应 2..6 体相互作用的耦合强度（可自由设定）
Kvec(2) = 0.5;    % 二体
Kvec(3) = 0;    % 三体
Kvec(4) = 25;    % 四体
% Kvec(5) = 0;    % 五体
% Kvec(6) = 0;    % 六体

dt = 0.1;         % 时间步长
T  = 20;           % 总时长
rng(2025);         % 频率扰动用的随机种子（可复现）

%% ========= 2) 提取网络数据 & 准备高阶度 =========
A = Sdata.A;              % n×n 邻接
S = Sdata.S;              % S{r}: r 个顶点的团列表（升序）
f = Sdata.f(:);
N = size(A,1);
Rmax = 4;  % 只用到 r<=6 的阶

% 计算（或读取）每个节点的各阶度 deg{r}(i)
deg = cell(Rmax,1);
if isfield(Sdata,'deg') && numel(Sdata.deg)>=Rmax
    for r = 2:Rmax, deg{r} = Sdata.deg{r}; end
else
    for r = 2:Rmax
        Er = S{r};
        if isempty(Er), deg{r} = zeros(N,1);
        else,           deg{r} = accumarray(Er(:), 1, [N,1], @sum, 0);
        end
    end
end

%% ========= 3) 构造自然频率 & 初相 =========
halfN  = floor(N/2);
omega  = zeros(N,1);
omega(1:halfN)        =  omega0 + sigma*rand(halfN,1);
omega(halfN+1:end)    = -omega0 - sigma*rand(N-halfN,1);

rng(2598);
theta1 = zeros(N/2,1)-2*pi/3*(rand(N/2,1)<0.1);  % 初相全 0
% unite1 = [0, 0, 0, 0, 0, 0, 0, 2*pi/3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];  % Fixed phase pattern with strict 9:1:0 ratio
% theta1 = repmat(unite1, 1, N/40);  % Initial phases for the first group
    % theta1 = (rand(N/2, 1) < 0.08) * 2/3*pi;  % (Alternative) Random initial condition
theta2 = flip(theta1);  % Initial phases for the second group (mirror of the first)
                         % The two groups have identical initial phases but asymmetric dynamics
theta0 = [theta1; theta2];  % Full phase vector

%% ========= 4) 时间推进（显式 Euler；如需更稳可换 ode45） =========
nT   = floor(T/dt) + 1;
t    = (0:nT-1).' * dt;
Theta = zeros(N, nT);
Theta(:,1) = theta0;
absz = zeros(1, nT);   % |z| 顺序参数模长
%% NEW: 两个半群的序参量（复数与模长）
z_half   = zeros(2, nT);         % z_half(1,:) 前半；z_half(2,:) 后半
absz_half= zeros(2, nT);         % 对应的模长 |z|

for k = 1:nT-1
    theta  = Theta(:,k);
    dtheta = omega;  % 频率基项

    % —— r 体耦合（r=2..6），仅沿网络中实际存在的 r 点团，并按节点 r 阶度归一化
    for r = 2:Rmax
        if Kvec(r)==0, continue; end
        Er = S{r};
        m  = size(Er,1);
        if m==0, continue; end

        sum_e    = sum(theta(Er), 2);             % 每条 r-超边的相位和（m×1）
        idx_all  = Er(:);                          % (m*r)×1：所有“超边-节点”展开
        sum_rep  = kron(ones(r,1), sum_e);        % (m*r)×1：重复的相位和
        theta_in = theta(idx_all);                % (m*r)×1：对应节点自身相位

        % 基本相互作用项：sin( sum_{j in e} θ_j - r * θ_i )
        contrib  = sin( sum_rep - r * theta_in );

        % —— 关键归一化：按“这个节点的 r 阶度”逐项除
        dr = deg{r}(idx_all);                     % (m*r)×1：每个项对应节点的 r 阶度
        % 对于存在于 r-超边中的节点，dr>=1；因此无需防 0
        contrib = contrib ./ dr;

        % 聚合到每个节点
        dtheta = dtheta + Kvec(r) * accumarray(idx_all, contrib, [N,1], @sum, 0);
    end

    % 一步 Euler
    Theta(:,k+1) = Theta(:,k) + dt * dtheta;

    % 顺序参数
    z = mean(exp(1i*Theta(:,k+1)));
    absz(1,k+1) = abs(z);
    %% NEW: 计算两半群的序参量
    idx1 = 1:halfN;                 % 前一半索引（与频率分配一致）
    idx2 = (halfN+1):N;             % 后一半索引
    z1 = mean(exp(1i*Theta(idx1, k+1)));
    z2 = mean(exp(1i*Theta(idx2, k+1)));
    z_half(1, k+1)    = z1;
    z_half(2, k+1)    = z2;
    absz_half(1, k+1) = abs(z1);
    absz_half(2, k+1) = abs(z2);

    
end

%% ========= 5) 简要结果 & 可视化 =========
% fprintf('Done. N=%d, steps=%d, f(0..%d)=%s\n', N, nT-1, numel(f)-1, mat2str(f.'));
% figure('Name','Order parameter |z|'); plot(t, absz, 'LineWidth',1.5); grid on
% xlabel('time'); ylabel('|z|'); title('|z| over time (higher-order Kuramoto with degree-normalization)');

% % 抽样画若干节点相位轨迹（包角到 [-pi,pi]，更易读）
% pick = min(25, N); idx = randperm(N, pick);
% figure('Name','Sample phases'); plot(t, wrapToPi(Theta(idx,:))'); grid on
% xlabel('time'); ylabel('\theta_i (wrapped)'); title(sprintf('Sample %d phases', pick));

%% NEW: 两个子图——前半/后半的序参量随时间
% —— 同一张图上用不同颜色绘制两半群的 |z|
figure('Name','Half-group order parameters (one axes)');
hold on; grid on
% 自定义两种容易区分的颜色（也可省略，使用 MATLAB 默认配色）
col1 = [0 0.4470 0.7410];   % 蓝
col2 = [0.8500 0.3250 0.0980]; % 橙

p1 = plot(t, absz_half(1,:), 'LineWidth', 1.8, 'Color', col1);
p2 = plot(t, absz_half(2,:), 'LineWidth', 1.8, 'Color', col2);

xlim([t(1) t(end)]);
ylim([0 1]);
xlabel('time'); ylabel('|z|');
title('前/后半振子序参量 |z|');
legend([p1 p2], {'前半','后半'}, 'Location','best');
hold off

%% ========= 辅助：自动找最近的 bundle =========
function latest = findLatestBundle()
    latest = "";
    d = dir('flag_n200_p0.500*');
    bestTime = -inf;
    for k = 1:numel(d)
        if ~d(k).isdir, continue; end
        m = dir(fullfile(d(k).name, 'flag_complex_bundle.mat'));
        if ~isempty(m)
            [~, idx] = max([m.datenum]);
            if m(idx).datenum > bestTime
                bestTime = m(idx).datenum;
                latest = fullfile(d(k).name, m(idx).name);
            end
        end
    end
end
