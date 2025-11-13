% ==========================================================
% scan_K1_parallel.m
% 在已构建的 flag 复形上做高阶 Kuramoto（至多 6 体）
% 并行扫描二体耦合 K1 ∈ [0,1]（常数，不爬坡），
% 记录并绘制末态两半群序参量模长之差。
% ==========================================================
clear; clc;

%% ========= 1) 选择要用的网络（.mat）与仿真参数 =========
bundleMat = "";  % 留空自动搜最新 flag_*/flagSF_*/flagWS_* 的 bundle
if bundleMat == "", bundleMat = string(findLatestBundle()); end
if bundleMat == "", error('未找到 flag_* / flagSF_* / flagWS_* 下的 flag_complex_bundle.mat'); end
Sdata = load(bundleMat);
fprintf('Loaded: %s\n', bundleMat);

% --- Kuramoto 参数（除 K1 外都固定）---
omega0 = 5;      % 基频幅值（前半 +omega0，后半 -omega0）
sigma  = 0.1;    % 频率扰动强度（与 rand 相乘）
Kvec_base = zeros(6,1);
Kvec_base(2) = 0;    % 二体耦合由扫描给出，这里占位
Kvec_base(3) = 0;    % 三体
Kvec_base(4) = 25;   % 四体
% Kvec_base(5) = 0;  % 五体（默认即 0）
% Kvec_base(6) = 0;  % 六体（默认即 0）
dt = 0.1;       % 时间步长
T  = 20;        % 总时长

% --- 扫描网格（可按需改精度）---
K1_list = linspace(0, 1, 21);  % 0:0.05:1
nK      = numel(K1_list);

%% ========= 2) 提取网络数据 & 准备高阶度 =========
A = Sdata.A;              %#ok<NASGU> % 如需可视化边层时可用
S = Sdata.S;              % S{r}: r 个顶点的团列表（升序）
f = Sdata.f(:);           %#ok<NASGU>
N = size(S{1},1);
Rmax = min(6, numel(S));  % 只用到网络实际提供的阶（且不超过 6）

% 读取或计算各阶“高阶度” deg{r}(i)
deg = cell(Rmax,1);
if isfield(Sdata,'deg') && numel(Sdata.deg) >= Rmax
    for r = 2:Rmax, deg{r} = Sdata.deg{r}; end
else
    for r = 2:Rmax
        Er = S{r};
        if isempty(Er), deg{r} = zeros(N,1);
        else,           deg{r} = accumarray(Er(:), 1, [N,1], @sum, 0);
        end
    end
end

%% ========= 3) 构造自然频率 & 初相（与之前一致）=========
halfN = floor(N/2);

rng(2025);                            % 频率扰动可复现
omega = zeros(N,1);
omega(1:halfN)       =  omega0 + sigma*rand(halfN,1);
omega(halfN+1:end)   = -omega0 - sigma*rand(N-halfN,1);

rng(2598);                            % 初相模式可复现
theta1 = zeros(halfN,1);
mask10 = (rand(halfN,1) < 0.10);      % 前半 10% 置为 2π/3，余为 0
theta1(mask10) = 2*pi/3;
theta2 = flip(theta1);                % 后半镜像
theta0 = [theta1; theta2];

idx_front = 1:halfN;
idx_back  = (halfN+1):N;

%% ========= 4) 并行扫描 K1（常数，不爬坡）=========
% 预分配输出
absz_front_end = zeros(1, nK);
absz_back_end  = zeros(1, nK);
Delta_end      = zeros(1, nK);

% 并行池
pool = gcp('nocreate');
if isempty(pool), parpool; end

nT = floor(T/dt) + 1;   % 时间步数（所有 K1 相同）

parfor kk = 1:nK
    K1 = K1_list(kk);
    % 固定其它阶耦合，替换二体为当前 K1
    Kvec = Kvec_base;
    Kvec(2) = K1;

    % 每个 K1 都从同一初相出发（不爬坡）
    theta = theta0;

    % 显式 Euler 推进
    for step = 1:nT-1
        dtheta = omega;  % 自然频率项

        % r 体相互作用：仅沿存在的 r 点团；按节点 r 阶度归一化
        for r = 2:Rmax
            Kr = Kvec(r);
            if Kr == 0, continue; end
            Er = S{r};
            m  = size(Er,1);
            if m == 0, continue; end

            sum_e    = sum(theta(Er), 2);             % 每条 r-超边相位和（m×1）
            idx_all  = Er(:);                          % (m*r)×1
            sum_rep  = kron(ones(r,1), sum_e);        % (m*r)×1
            theta_in = theta(idx_all);                % (m*r)×1

            contrib  = sin( sum_rep - r * theta_in ); % 基本相互作用
            dr       = deg{r}(idx_all);               % 归一化：节点的 r 阶度
            contrib  = contrib ./ dr;

            dtheta = dtheta + Kr * accumarray(idx_all, contrib, [N,1], @sum, 0);
        end

        theta = theta + dt * dtheta;
    end

    % 末态两半群序参量
    z_front = mean( exp(1i * theta(idx_front)) );
    z_back  = mean( exp(1i * theta(idx_back)) );

    absz_front_end(kk) = abs(z_front);
    absz_back_end(kk)  = abs(z_back);
    Delta_end(kk)      = abs( absz_front_end(kk) - absz_back_end(kk) );
end

% %% ========= 5) 作图：末态两半群 |z| 之差 vs K1 =========
% figure('Name','End-state |z_front|-|z_back| vs K1');
% plot(K1_list, Delta_end, '-o', 'LineWidth', 1.8); grid on
% xlabel('K_1 (two-body coupling)');
% ylabel('||z_{front}| - |z_{back}|| (end state)');
% title('末态两半群序参量模长之差 vs 二体耦合 K_1');
% 
% % 可选：同时画两半群的末态 |z|
% figure('Name','End-state |z| (front/back) vs K1');
% plot(K1_list, absz_front_end, '-o', 'LineWidth', 1.6); hold on;
% plot(K1_list, absz_back_end,  '-s', 'LineWidth', 1.6); grid on; hold off;
% xlabel('K_1'); ylabel('|z| (end)');
% legend('前半','后半','Location','best');
% title('末态两半群 |z| vs K_1');

%% ========= 6) （可选）保存结果 =========
out = struct();
out.bundleMat = bundleMat;
out.K1_list   = K1_list;
out.absz_front_end = absz_front_end;
out.absz_back_end  = absz_back_end;
out.Delta_end      = Delta_end;
save('scan_K1_parallel_results.mat','-struct','out');

fprintf('Done. 结果已保存到 scan_K1_parallel_results.mat\n');

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
