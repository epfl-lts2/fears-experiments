%VARIOUS_GRAPHS_PROJ Function to perform eigenapproximation with graphs

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

function [ data ] = various_graphs_proj( params )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    data = {};
    data.params = [];
    data.measures = [];
    
    params.k = params.k_large;
    
    %Sensor network
    paramsensor.connected = 1;
    g11 = @(N)gsp_sensor(params.N, paramsensor);
    measures = compute_proj(g11, params);

    params.name = 'sensor';
    data.params = [data.params, params];
    data.measures = [data.measures, measures];

    %SBM
    sbmparams.p = 16*params.k / params.N;
    sbmparams.q = 4 / params.N;
    g12 = @(N)gsp_stochastic_block_graph(N, params.k, sbmparams);
    measures = compute_proj(g12, params);

    params.name = 'sbm';
    data.params = [data.params, params];
    data.measures = [data.measures, measures];
    
    %Minnesota road network
    g14 = @(N)gsp_minnesota();
    measures = compute_proj(g14, params);

    params.name = 'minnesota';
    data.params = [data.params, params];
    data.measures = [data.measures, measures];

    params.k = params.k_small;
    
    %Swiss roll
    g21 = @(N)gsp_swiss_roll(N,randi(10000));
    measures = compute_proj(g21, params);

    params.name = 'swiss';
    data.params = [data.params, params];
    data.measures = [data.measures, measures];

    %Stanford Bunny
    g22 = @(N)gsp_bunny;
    tmpN = params.N;
    params.N = 2500;
    measures = compute_proj(g22, params);

    params.name = 'bunny';
    data.params = [data.params, params];
    params.N = tmpN;
    data.measures = [data.measures, measures];

    %Image patch graph
    img = imread('data/barbara.png');
    img = imresize(img, 0.25);
    parampatch.rho = 50;
    g33 = @(N)gsp_patch_graph(img, parampatch);
    tmpN = params.N;
    params.N = 16384;
    measures = compute_proj(g33, params);

    params.name = 'image';

    data.params = [data.params, params];
    params.N = tmpN;
    data.measures = [data.measures, measures];

end

