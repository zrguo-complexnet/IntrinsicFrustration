% ==========================================================
% build_flag_complex_and_save.m
% 先构建随机团（flag）复形，导出各阶连边数据
% 生成：.mat 总包 + 各阶 CSV + 简要 summary.txt
% ==========================================================
clear; clc;

%% ---- 参数区（按需修改）----
n       = 200;        % 节点数
p       = 0.5;       % 边出现概率（G(n,p)）
maxDim  = 4;          % 最高维：4 表示最多找 5 个点的团
rngSeed = 2025;       % 随机种子，保证可重复

% 可选：是否构建“节点-超边”关联矩阵 H_r（稀疏存 mat；CSV 不适合）
% 出于内存考虑，默认只到 r=4（即 4 个点的团=三维单纯形）
buildIncidenceUpTo = 4;

%% ---- 1) 构建随机团（flag）复形 ----
% 需要 randomFlagComplex.m（已给出）
[A, S, f] = randomFlagComplex(n, p, maxDim, rngSeed);

% 实际最高 r（可能没有高维团，S 会提前结束）
Rmax = numel(S);

% 每阶“节点入度”：节点属于多少个 r 点团
deg = cell(Rmax,1);
for r = 2:Rmax
    Er = S{r};
    if ~isempty(Er)
        deg{r} = accumarray(Er(:), 1, [n,1], @sum, 0);
    else
        deg{r} = zeros(n,1);
    end
end

%% ---- 2) 构建（可选）节点-超边关联矩阵 H_r（稀疏，便于后续计算）----
H = cell(Rmax,1);  % H{r} 尺寸：N × (#r点团)，条目为 0/1（节点是否属于该超边）
for r = 2:min(Rmax, buildIncidenceUpTo)
    Er = S{r};
    m  = size(Er,1);
    if m == 0, H{r} = sparse(n,0); continue; end
    rows = reshape(Er.', [], 1);            % r*m × 1（每条超边的 r 个节点）
    cols = repelem((1:m).', r, 1);          % r*m × 1（对应的超边列号）
    vals = ones(numel(rows),1);
    H{r} = sparse(rows, cols, vals, n, m);  % 稀疏关联矩阵
end

%% ---- 3) 组织输出目录+文件名 ----
stamp  = datestr(now, 'yyyymmdd_HHMMSS');
outdir = sprintf('flag_n%d_p%.3f_d%d_%s', n, p, maxDim, stamp);
if ~exist(outdir, 'dir'), mkdir(outdir); end

matfile = fullfile(outdir, 'flag_complex_bundle.mat');

%% ---- 4) 保存 .mat 总包（建议后续加载此文件进行模拟）----
meta = struct('n',n,'p',p,'maxDim',maxDim,'rngSeed',rngSeed,...
              'createdAt', datestr(now,'yyyy-mm-ddTHH:MM:SS'),...
              'note','Indices are 1-based (MATLAB). S{r} stores r-vertex cliques.');
save(matfile, 'A','S','f','deg','H','meta','-v7.3');  % -v7.3 适合大稀疏矩阵

%% ---- 5) 保存各阶 CSV（每行一条超边/团；列是 r 个节点编号，升序）----
% 顶点（S{1}）通常无需导出为 CSV；从 r=2 开始
for r = 2:Rmax
    Er = S{r};
    if isempty(Er), continue; end
    csv_r = fullfile(outdir, sprintf('simplices_r%d.csv', r));
    writematrix(Er, csv_r);   % 默认逗号分隔
end

% 同时导出各阶节点入度（每个节点在该阶出现的次数）
for r = 2:Rmax
    dr = deg{r};
    if isempty(dr), continue; end
    csv_d = fullfile(outdir, sprintf('degree_r%d.csv', r));
    writematrix(dr, csv_d);
end

% 边的邻接矩阵也导出一份（0/1），便于外部工具使用
csv_A = fullfile(outdir, 'adjacency_edges_r2.csv');
writematrix(double(A), csv_A);

%% ---- 6) 写一个简要 summary.txt（方便快速查看规模）----
summaryFile = fullfile(outdir, 'summary.txt');
fid = fopen(summaryFile, 'w');
fprintf(fid, 'Random flag complex summary\n');
fprintf(fid, 'n=%d, p=%.6f, maxDim=%d, seed=%d\n', n,p,maxDim,rngSeed);
fprintf(fid, 'f-vector (0..%d dim): ', Rmax-1);
fprintf(fid, '%s\n', mat2str(f(:).'));
for r = 2:Rmax
    fprintf(fid, 'r=%d  (#%d-vertex cliques) = %d\n', r, r, size(S{r},1));
end
fprintf(fid, '\nFiles:\n');
fprintf(fid, '  %s\n', matfile);
fprintf(fid, '  %s\n', csv_A);
for r = 2:Rmax
    if ~isempty(S{r})
        fprintf(fid, '  simplices_r%d.csv\n', r);
        fprintf(fid, '  degree_r%d.csv\n',   r);
    end
end
fclose(fid);

%% ---- 7) 终端提示 ----
fprintf('Done. Output directory: %s\n', outdir);
fprintf('Summary: %s\n', summaryFile);
