%% ============================================================
%  Main script: Load structural connectivity, extract 4-cliques,
%  build node–to–clique incidence table, and run high-order 
%  Kuramoto dynamics (pairwise + 4-body coupling).
%
%  This script assumes that CIJ_resampled_average is already 
%  prepared and stored in a .mat file.
%% ============================================================

%% 0. Clear workspace and load the structural connectivity matrix
clear; clc;

load('your_data.mat', 'CIJ_resampled_average');   % Replace with your file name

% Ensure a binary, symmetric adjacency matrix with zero diagonal
A = (CIJ_resampled_average ~= 0);
A = triu(A, 1);
A = A + A.';
A(1:size(A,1)+1:end) = 0;

N = size(A, 1);   % Number of nodes

%% 1. Extract 4-cliques from the adjacency matrix
%    (find_cliques returns both 3- and 4-cliques; only the latter is needed)
[~, four_cliques] = find_cliques(A);

if isempty(four_cliques)
    error('No 4-cliques detected in the given network.');
end

%% 2. Construct node–to–clique incidence matrix
%    clique_matrix(i,:) lists all 4-cliques that node i belongs to;
%    max_cliques(i) is the number of such cliques.
[clique_matrix, max_cliques] = organize_4cliques(four_cliques);

% Sanity check: clique_matrix must have one row per node
if size(clique_matrix, 1) ~= N
    error('Dimension mismatch: clique_matrix has %d rows, CIJ has %d nodes.', ...
          size(clique_matrix,1), N);
end

%% 3. Set Kuramoto parameters and run the dynamics
K2 = 15;            % 4-body coupling strength
K1 = 3;             % Pairwise coupling strength
num_steps = 5000;   % Total simulation steps
dt = 0.01;          % Time step

kuramoto_simulation( ...
        clique_matrix, ...      % Node–to–4-clique table
        four_cliques, ...       % Raw 4-clique list
        max_cliques, ...        % Degree in 4-clique layer
        K1, K2, ...             % Pairwise and 4-body couplings
        num_steps, dt, ...      % Simulation parameters
        CIJ_resampled_average); % Pairwise adjacency matrix
