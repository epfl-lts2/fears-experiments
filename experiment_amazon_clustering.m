%EXPERIMENT_AMAZON_CLUSTERING Script for amazon clustering

% Copyright (C) 2016 Johan Paratte, Lionel Martin.
% This file is part of the Reproducible Research code base for the Fast
% Eigenspace Approximation using Random Signals (FEARS) method.
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% If you use this code please kindly cite
%     Paratte, Johan, and Lionel Martin. 
%     "Fast Eigenspace Approximation using Random Signals."
%     arXiv preprint arXiv:1611.00938 (2016).
% https://arxiv.org/abs/1611.00938

load 'data/com-amazon.ungraph2.txt';
N = max(com_amazon_ungraph2(:, 2));
S = sparse(com_amazon_ungraph2(:, 1), com_amazon_ungraph2(:, 2), ones(size(com_amazon_ungraph2, 1), 1), N, N);
W = S + S';
null_idx = find(sum(W) == 0);
W(null_idx, :) = [];
W(:, null_idx) = [];

G = gsp_graph(W);
G = gsp_create_laplacian(G, 'normalized');
G = gsp_estimate_lmax(G);

k_values = [250, 500, 1000];
order = 50;

nb_exp = numel(k_values); 

modularities = zeros(nb_exp, 4);
times = zeros(nb_exp, 4, 2);
measures = zeros(nb_exp, 3, 2);

for idx_k=1:nb_exp
    k = k_values(idx_k); G.k = k;
    fprintf('Clustering of k = %d clusters\n', k);
    params.order = order; params.verbose = 1;
    
    % 1. Us - kmeans
    fprintf('Starting our filtering...\n'); tic;
    [basis, ~, out_params] = gsp_eigenspace_estimation(G, k, params);
    T_filt = toc;
    times(idx_k, 1, 1) = T_filt;
    times(idx_k, 2, 1) = T_filt;
    fprintf('Our filtering completed (time: %f).\nStarting kmeans with these features... ', T_filt);

    tic;
    assign_us = kmeans(basis, k, 'Replicates', 3, 'MaxIter', 1000, 'Options', statset('UseParallel',1));
    times(idx_k, 1, 2) = toc;

    modularities(idx_k, 1) = compute_modularity(assign_us, G.W);
    fprintf('Our method terminated (time: %f). Modularity: %f.\nStarting compressive kmeans on our features...\n', toc, modularities(idx_k, 1));


    % 2. Us - compressive
    tic;
    assign_CUS = compressive_kmeans(G, basis, out_params.lk, order);
    times(idx_k, 2, 2) = toc;

    modularities(idx_k, 2) = compute_modularity(assign_CUS, G.W);
    fprintf('Our method terminated (time: %f). Modularity: %f.\nStarting eigs... ', toc, modularities(idx_k, 2));

    % 3. SC
    tic;
    [Uk, ek] = eigs(G.L, k, 'sm');
    Uk = fliplr(Uk); ek = flip(diag(ek));
    times(idx_k, 3, 1) = toc;
    fprintf('Eigs completed (time: %f).\nStarting kmeans with these features (SC)... ', toc);

    tic;
    assign_SC = kmeans(Uk, k, 'Replicates', 3, 'MaxIter', 1000, 'Options', statset('UseParallel',1));
    times(idx_k, 3, 2) = toc;
    
    modularities(idx_k, 3) = compute_modularity(assign_SC, G.W);
    fprintf('SC terminated (time: %f). Modularity: %f.\nStarting CSC features extraction... ', toc, modularities(idx_k, 3));

    % 4. CSC
    tic;
    param_CSC.poly_order = order; param_CSC.only_features = 1;
    [CSC_feat, CSC_lk, ~, ~, ~, ~] = CSC(G, param_CSC);
    times(idx_k, 4, 1) = toc;
    fprintf('CSC features extracted (time: %f).\nStarting compressive kmeans with these features (CSC)... ', toc);

    tic;
    assign_CSC = compressive_kmeans(G, CSC_feat, CSC_lk, order);
    times(idx_k, 4, 2) = toc;
    
    modularities(idx_k, 4) = compute_modularity(assign_CSC, G.W);
    fprintf('CSC terminated (time: %f). Modularity: %f.\n\n', toc, modularities(idx_k, 4));
    
    % Measures
    methods = [assign_us, assign_CUS, assign_CSC];
    parfor id = 1:3
        measures(idx_k, id, 1) = rand_index(assign_SC, methods(:, id));
    end

    parfor id = 1:3
        measures(idx_k, id, 2) = compute_nmi(assign_SC, methods(:, id));
    end
end