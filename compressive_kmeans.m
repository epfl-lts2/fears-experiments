function [assign] = compressive_kmeans(G, features, lk, order)
%COMPRESSIVE_KMEANS Function which computes compressive kmeans

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

    weight = ones(G.N, 1)/G.N;
    n = round(4 * G.k * log(G.k));
    ind_obs = datasample(1:G.N, n, 'Replace', false, 'Weights', weight);
    X_lk_est_DS = features(ind_obs, :);
    IDX_LD = kmeans(X_lk_est_DS , G.k, 'Replicates', 20, 'Options', statset('UseParallel',1));
    C_obs_LD = sparse(1:n, IDX_LD, 1, n, G.k);
    [~,JCH_HP] = jackson_cheby_poly_coefficients(lk, G.lmax, [0, G.lmax], order);
    C_est = zeros(G.N, size(C_obs_LD,2));
    parfor ki=1:size(C_obs_LD,2)
        c_obs = C_obs_LD(:,ki);
        C_est(:, ki) = interpolate_on_complete_graph(c_obs, ind_obs, @(x)gsp_cheby_op(G, JCH_HP, x), 1e-3, G.N, 'gmres');
    end
    [~, assign] = max(C_est./repmat(sqrt(sum(C_est.^2, 1)), G.N, 1), [], 2);
end