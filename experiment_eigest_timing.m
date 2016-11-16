%EXPERIMENT_EIGEST_TIMING Script for timing the eigenestimation experiments

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


%% Experiment 1

N = 10000;

G = gsp_sensor(N);
G = gsp_estimate_lmax(G);

%ks = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048];

ks = [2, 4, 8, 16, 32, 64];

timing = zeros(length(ks), 4);

large_scale = 0;

for kk = 1:length(ks)
    data = timing_all(G, ks(kk), large_scale);
    timing(kk, 1) = data.eigs;
    timing(kk, 2) = data.our;
    timing(kk, 3) = data.power;
    timing(kk, 4) = data.csc;
end

figure;
plot(ks, timing(:,1))
hold on; plot(ks, timing(:,2))
hold on; plot(ks, timing(:,3))
hold on; plot(ks, timing(:,4))
legend('eigs', 'our', 'power', 'csc');

%% Experiment 2

N = 10000;

G = gsp_sensor(N);
G = gsp_estimate_lmax(G);

ks = [32, 64, 128, 256, 512, 1024, 2048];

timing = zeros(length(ks), 3);

large_scale = 1;

for kk = 1:length(ks)
    data = timing_all(G, ks(kk), large_scale);
    timing(kk, 1) = data.eigs;
    timing(kk, 2) = data.our;
    timing(kk, 3) = data.csc;
end

figure;
plot(ks, timing(:,1))
hold on; plot(ks, timing(:,2))
hold on; plot(ks, timing(:,3))
legend('eigs', 'our', 'csc');

%% Expriment 3 with k fixed, N moving

k = 250;
ks = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048];
ns = ks * 1000;

large_scale = 1;

timing = zeros(length(ns), 3);

for nn = 1:length(ns)
    n = ns(nn);
    clear G;
    param.sensor.k = 10;
    G = gsp_sensor(n, param.sensor);
    G = gsp_estimate_lmax(G);
    
    data = timing_all(G, k, large_scale);
    timing(nn, 1) = data.eigs;
    timing(nn, 2) = data.our;
    timing(nn, 3) = data.csc;
end


%% Expriment 3  with k log(N), N moving

%k = 250;
ks = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048];
ns = ks * 1000;

large_scale = 1;

timing = zeros(length(ns), 3);

for nn = 1:length(ns)
    
    n = ns(nn);
    k = round(log(n));
    clear G;
    param.sensor.k = 10;
    G = gsp_sensor(n, param.sensor);
    G = gsp_estimate_lmax(G);
    
    data = timing_all(G, k, large_scale);
    timing(nn, 1) = data.eigs;
    timing(nn, 2) = data.our;
    timing(nn, 3) = data.csc;
    
end

%% Expriment 4  with k sqrt(N), N moving

%k = 250;
ks = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048];
ns = ks * 1000;

large_scale = 1;

timing = zeros(length(ns), 3);

for nn = 1:length(ns)
    
    n = ns(nn);
    k = round(0.5*sqrt(n));
    clear G;
    param.sensor.k = 10;
    G = gsp_sensor(n, param.sensor);
    G = gsp_estimate_lmax(G);
    
    data = timing_all(G, k, large_scale);
    timing(nn, 1) = data.eigs;
    timing(nn, 2) = data.our;
    timing(nn, 3) = data.cscs;
end

%%
fig = gcf;
cmaps = cbrSelector('name','Set1', 'length', 8, 'plot', false);
cmap = uint8(cmaps{1});
ax = fig.CurrentAxes;
plot(ns, timing(:,1), 'color', cmap(1,:));
hold on; plot(ns, timing(:,2), 'color', cmap(2,:));
hold on; plot(ns, timing(:,3), 'color', cmap(3,:));
legend('eigs', 'our', 'csc');
