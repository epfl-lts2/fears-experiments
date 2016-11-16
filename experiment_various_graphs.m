%EXPERIMENT_VARIOUS_GRAPHS Script for FEARS evaluation on various graphs

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

%% Quality of approximation with respect various graphs

params.N = 10000;
params.k_large = 25;
params.k_small = 25;
params.nb_avg = 50;
params.order = 500;
%params.verbose = 1;

% Do the experiment
ts = tic;
data2 = various_graphs_proj(params);
te = toc(ts);

% Write results
print_proj_results(data2);