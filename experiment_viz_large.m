%EXPERIMENT_VIZ_LARGE Script used to perform largescale visualization

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

%% Init

gsp_start;

%% Import data
% Note : LiveJournal is to big to be included in the data available here
load 'data/mnist.mat';

paramnn.k = 10;
paramnn.use_flann = 1;

G = gsp_nn_graph(data', paramnn);
G = gsp_estimate_lmax(G);
G = gsp_create_laplacian(G, 'normalized');

mksz = 8;
k = 2;
l1 = labels;

%% Laplacian eigenmaps

t1 = tic;
[U, V] = eigs(G.L, diag(G.d), k+1, 'sm');
t2 = toc(t1);
fprintf('Time for laplacian eigenmaps %f s\n', t2);

figure;
subplot(221)
scatter(U(:,1), U(:,2), mksz, l1);
title('Laplacian eigenmaps');
colormap(jet)

%% FEARS

params = {};
params.order = 250;

t1 = tic;
[Bk, approx_U] = gsp_eigenspace_estimation(G, k, params);

t2 = toc(t1);
fprintf('Time for our method %f s\n', t2);

% figure;
% scatter(approx_U(:,1), approx_U(:,2), mksz, l1');
% title('Our method M');
% colormap(jet)

figure;
subplot(222)
scatter(Bk(:,1), Bk(:,2), mksz, l1');
title('Our method Bk');
colormap(jet)

%% T-SNE
% You need to download and compile BH t-SNE to run this experiment
% available here : https://github.com/lvdmaaten/bhtsne

t1 = tic;
mappedX = fast_tsne(data', k);
t2 = toc(t1);
fprintf('Time for t-SNE %f s\n', t2);

%figure;
subplot(223);
scatter(mappedX(:,1), mappedX(:,2), mksz, l1');
title('t-SNE');
colormap(jet)
