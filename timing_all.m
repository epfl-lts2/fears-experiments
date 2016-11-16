function [ data ] = timing_all(G, k, fast)
%TIMING_ALL Function which compute timing for different methods

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

    %Eigs
    t1 = tic;
    [V, D] = eigs(G.L, k, 'sm');
    tim = toc(t1);
    data.eigs = tim;
    
    %FEARS
    clear params_filt;
    params_filt.is_exact = 0; % use exact lowpass filtering or polynomial approximation
    params_filt.order = 50; % order of the polynomial approx
    t1 = tic;
    [Bk, approx_U, ~] = gsp_eigenspace_estimation(G, k, params_filt);
    tim = toc(t1);
    
    data.our = tim;
    
    %Power method
    if ~fast
        t1 = tic;
        Wk = (speye(G.N) - G.L);
        approx_U = (Wk^k)*randn(G.N, k);
        [Bk, ~, ~] = svd(approx_U, 'econ');
        tim = toc(t1);
        data.power = tim;
    end
    
    % CSC
    param_CSC.poly_order = params_filt.order;
    param_CSC.only_features = 1;
    param_CSC.lap_type = 'normalized';
    G.k = k;
    t1 = tic;
    CSC(G,param_CSC);
    tim = toc(t1);
    
    data.csc = tim;
end

