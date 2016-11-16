function [  ] = print_proj_results( data )
%PRINT_PROJ_RESULTS Function which formats the data from various graphs

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

    fprintf('=======================================================\n');
    fprintf('                        Results                        \n');
    fprintf('=======================================================\n');

    N = length(data.params);
    fprintf('name | ');
    for ii = 1:N
       fprintf('%s \t', data.params(ii).name); 
    end
    fprintf('\n');
    fprintf('-------------------------------------------------------\n');

    fprintf('N    | ');
    for ii = 1:N
       fprintf('%d \t', data.params(ii).N); 
    end
    fprintf('\n');

    fprintf('k    | ');
    for ii = 1:N
       fprintf('%d \t', data.params(ii).k); 
    end
    fprintf('\n');

    fprintf('mp   | ');
    tot = 0;
    for ii = 1:N
       fprintf('%0.2f $\\pm$ %0.2f \t', data.measures(ii).mp, data.measures(ii).mpstd);
       tot = tot + data.measures(ii).mp;
    end
    fprintf('| %0.2f', tot/N);
    fprintf('\n');

    fprintf('mp-s | ');
    tot = 0;
    for ii = 1:N
       fprintf('%0.2f $\\pm$ %0.2f \t', data.measures(ii).mps, data.measures(ii).mpsstd);
       tot = tot + data.measures(ii).mps;
    end
    fprintf('| %0.2f', tot / N);
    fprintf('\n');

    fprintf('mp-f | ');
    tot = 0;
    for ii = 1:N
       fprintf('%0.2f $\\pm$ %0.2f \t', data.measures(ii).mpf, data.measures(ii).mpfstd);
        tot = tot + data.measures(ii).mpf;
    end
    fprintf('| %0.2f', tot / N);
    fprintf('\n');

    fprintf('it-s | ');
    tot = 0;
    for ii = 1:N
       fprintf('%0.2f $\\pm$ %0.2f \t', data.measures(ii).mis, data.measures(ii).misstd);
       tot = tot + data.measures(ii).mis;
    end
    fprintf('| %0.2f', tot / N);
    fprintf('\n');


    fprintf('it-f | ');
    tot = 0;
    for ii = 1:N
        fprintf('%0.2f $\\pm$ %0.2f \t', data.measures(ii).mif, data.measures(ii).mifstd);
        tot = tot + data.measures(ii).mif;
    end
    fprintf('| %0.2f', tot/N);
    fprintf('\n');

    fprintf('kd-s | ');
    tot = 0;
    for ii = 1:N
       fprintf('%0.2f $\\pm$ %0.2f \t', data.measures(ii).mks, data.measures(ii).mksstd);
       tot = tot + data.measures(ii).mks;
    end
    fprintf('| %0.2f', tot/N);
    fprintf('\n');


    fprintf('kd-f | ');
    tot = 0;
    for ii = 1:N
       fprintf('%0.2f $\\pm$ %0.2f \t', data.measures(ii).mkf, data.measures(ii).mkfstd);
        tot = tot + data.measures(ii).mkf;
    end
    fprintf('| %0.2f', tot / N);
    fprintf('\n');
    
    fprintf('lkd-s | ');
    tot = 0;
    for ii = 1:N
       fprintf('%0.2f $\\pm$ %0.2f \t', data.measures(ii).lks, data.measures(ii).lksstd);
       tot = tot + data.measures(ii).lks;
    end
    fprintf('| %0.2f', tot/N);
    fprintf('\n');


    fprintf('lkd-f | ');
    tot = 0;
    for ii = 1:N
        fprintf('%0.2f $\\pm$ %0.2f \t', data.measures(ii).lkf, data.measures(ii).lkfstd);
        tot = tot + data.measures(ii).lkf;
    end
    fprintf('| %0.2f', tot / N);
    fprintf('\n');


    fprintf('=======================================================\n');

end

