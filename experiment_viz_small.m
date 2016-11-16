%EXPERIMENT_VIZ_SMALL Script used to perform smallscale visualization

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

%% Data init
clear G;

N = 10000;
params_sr.k = 10;
[G, tt] = gsp_swiss_roll(N, 0, params_sr);
G = gsp_estimate_lmax(G);
G = gsp_create_laplacian(G, 'normalized');

k = 2;
t1 = tt;

gsp_plot_signal(G, t1);

%% Laplacian eigenmaps

mksz = 8;

[U, V] = eigs(G.L, diag(G.d), k+1, 'sm');

scatter(U(:,1), U(:,2), mksz, t1');
colormap(jet)

%% FEARS

params = {};
params.order = 1000;
params.verbose = 1;
params.fast_lk = 1;

G = gsp_estimate_lmax(G);
G.lmax = G.lmax / 1.01;

[Bk, approx_U] = gsp_eigenspace_estimation(G, k, params);

subplot(223);
scatter(approx_U(:,1), approx_U(:,2), mksz, t1');
title('Our method M');
colormap(jet)

subplot(224);
scatter(Bk(:,1), Bk(:,2), mksz, t1');
title('Our method Bk');
colormap(jet)
