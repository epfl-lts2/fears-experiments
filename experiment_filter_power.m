%EXPERIMENT_FILTER_POWER Script for the experiment on filter power

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

%% Study the effect of the polynomial power on the approximation

N = 500;

G = gsp_sensor(N);
G = gsp_estimate_lmax(G);
G = gsp_compute_fourier_basis(G);
Uk = G.U(:,1:k);


%% Cheby vs Jackson-Cheby
m = 300;
k = 90;
lk = G.e(k);


%figure;
paramplot.plot_eigenvalues = 1;
paramplot.cla = 0;
h = @(x) x < lk;
gsp_plot_filter(G, h, paramplot);


paramplot.plot_eigenvalues = 0; %avoid redrawing eigenvalues
cheby_approx = gsp_approx_filter(G, h, m);
hold on;
gsp_plot_filter(G, cheby_approx, paramplot);

jch_approx = approx_filter_jch(G, lk, m, N);
hold on;
gsp_plot_filter(G, jch_approx, paramplot);



%% Jackson-Cheby increasing m
k = 90;
lk = G.e(k);

ms = [50, 100, 200, 400];

%figure;
paramplot.plot_eigenvalues = 1;
paramplot.cla = 0;

h = @(x) x < lk;
gsp_plot_filter(G, h, paramplot);

for m = ms

    jch_approx = approx_filter_jch(G, lk, m, N);
    hold on;
    gsp_plot_filter(G, jch_approx, paramplot);

end