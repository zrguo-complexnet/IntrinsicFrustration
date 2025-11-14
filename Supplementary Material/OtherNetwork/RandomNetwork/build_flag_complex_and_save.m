% ==========================================================
% build_flag_complex_and_save.m
% Constructs a random flag complex G(n,p) and exports:
%   - All simplices S{r}, adjacency matrix A
%   - Optional incidence matrices H_r (sparse)
%   - Per-dimension CSV files and a summary report
%   - A consolidated .mat bundle for later simulations
% ==========================================================
clear; clc;

%% ---------------------- Parameter settings ----------------------
n       = 200;        % Number of vertices
p       = 0.5;        % Edge probability for G(n,p)
maxDim  = 4;          % Maximum simplex dimension (r-vertex cliques with r ≤ maxDim+1)
rngSeed = 2025;       % Random seed for reproducibility

% Optional: construct node–hyperedge incidence matrices H_r (sparse)
% For memory efficiency, typically only up to r = 4
buildIncidenceUpTo = 4;

%% ---------------------- 1) Construct flag complex ----------------------
% Requires the function randomFlagComplex.m
[A, S, f] = randomFlagComplex(n, p, maxDim, rngSeed);

% Actual highest r (may be smaller if no higher-dimensional cliques exist)
Rmax = numel(S);

% Node degrees per dimension: number of r-vertex cliques containing each node
deg = cell(Rmax,1);
for r = 2:Rmax
    Er = S{r};
    if ~isempty(Er)
        deg{r} = accumarray(Er(:), 1, [n,1], @sum, 0);
    else
        deg{r} = zeros(n,1);
    end
end

%% ---------------------- 2) Optional incidence matrices ----------------------
% H{r} is an n × (#r-simplices) sparse matrix with entries 0/1
H = cell(Rmax,1);
for r = 2:min(Rmax, buildIncidenceUpTo)
    Er = S{r};
    m  = size(Er,1);
    if m == 0
        H{r} = sparse(n,0);
        continue;
    end
    rows = reshape(Er.', [], 1);
    cols = repelem((1:m).', r, 1);
    vals = ones(numel(rows),1);
    H{r} = sparse(rows, cols, vals, n, m);
end

%% ---------------------- 3) Output directory ----------------------
stamp  = datestr(now, 'yyyymmdd_HHMMSS');
outdir = sprintf('flag_n%d_p%.3f_d%d_%s', n, p, maxDim, stamp);
if ~exist(outdir, 'dir'), mkdir(outdir); end

matfile = fullfile(outdir, 'flag_complex_bundle.mat');

%% ---------------------- 4) Save consolidated .mat bundle ----------------------
meta = struct('n',n,'p',p,'maxDim',maxDim,'rngSeed',rngSeed, ...
              'createdAt', datestr(now,'yyyy-mm-ddTHH:MM:SS'), ...
              'note','Indices are 1-based. S{r} stores all r-vertex cliques.');
save(matfile, 'A','S','f','deg','H','meta','-v7.3');

%% ---------------------- 5) Export per-dimension CSV files ----------------------
% Export r-vertex cliques (S{1} skipped intentionally)
for r = 2:Rmax
    Er = S{r};
    if isempty(Er), continue; end
    csv_r = fullfile(outdir, sprintf('simplices_r%d.csv', r));
    writematrix(Er, csv_r);
end

% Export node participation counts in each dimension
for r = 2:Rmax
    dr = deg{r};
    if isempty(dr), continue; end
    csv_d = fullfile(outdir, sprintf('degree_r%d.csv', r));
    writematrix(dr, csv_d);
end

% Export adjacency matrix (0/1)
csv_A = fullfile(outdir, 'adjacency_edges_r2.csv');
writematrix(double(A), csv_A);

%% ---------------------- 6) Write summary report ----------------------
summaryFile = fullfile(outdir, 'summary.txt');
fid = fopen(summaryFile, 'w');
fprintf(fid, 'Random Flag Complex Summary\n');
fprintf(fid, 'n=%d, p=%.6f, maxDim=%d, seed=%d\n', n,p,maxDim,rngSeed);
fprintf(fid, 'f-vector (0..%d dim): %s\n', Rmax-1, mat2str(f(:).'));
for r = 2:Rmax
    fprintf(fid, 'r=%d:  #(%d-vertex cliques) = %d\n', r, r, size(S{r},1));
end
fprintf(fid, '\nFiles:\n  %s\n  %s\n', matfile, csv_A);
for r = 2:Rmax
    if ~isempty(S{r})
        fprintf(fid, '  simplices_r%d.csv\n', r);
        fprintf(fid, '  degree_r%d.csv\n',   r);
    end
end
fclose(fid);

%% ---------------------- 7) Terminal message ----------------------
fprintf('Done. Output directory: %s\n', outdir);
fprintf('Summary saved to: %s\n', summaryFile);
