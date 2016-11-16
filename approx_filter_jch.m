function gt = approx_filter_jch(G, lk, m, N, param)
%APPROX_FILTER_JCH Function which gets the approximation of a jch poly
%inspired by gsp_approx_filter

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

if nargin < 5
    param = struct;
end

if ~isfield(param,'verbose'), param.verbose = 1; end;

if nargin < 3
    m = 30;
end

if nargin < 4
   N = m+1; 
end

if isstruct(G)
    if ~isfield(G,'lmax');
        G = gsp_estimate_lmax(G);
        if param.verbose
        warning(['GSP_APPROX_FILTER: The variable lmax is not ',...
            'available. The function will compute it for you. ',...
            'However, if you apply many time this function, you ',...
            'should precompute it using the function: ',...
            'gsp_estimate_lmax']);
        end
    end
  arange = [0, G.lmax];
else
  arange = G;
end
  

[~, jch] = jackson_cheby_poly_coefficients(0, lk, [0, G.lmax], m);
gt = @(x) gsp_cheby_eval(x,jch,arange);

% Nf = length(g);
% gt = cell(Nf,1);
% 
% for ii = 1:Nf
%     gt{ii} = @(x) gsp_cheby_eval(x,jch,arange);
% end
% 
% if Nf == 1;
%     gt = gt{1};
% end



end
