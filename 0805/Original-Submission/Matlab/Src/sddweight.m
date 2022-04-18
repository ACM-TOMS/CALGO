function [dsav, xsav, ysav, itssav, rhosav, iitssav] = ...
    sddweight(A, W, kmax, alphamin, lmax, rhomin, yinit)

%SDDWEIGHT  Weighted Semidiscrete Decomposition.
%
%   [D, X, Y] = SDDWEIGHT(A, W) produces discrete matrices X and Y and a
%   vector D such that X *diag(D) * Y' is the weighted SDD of rank 100 that
%   approximates A in the norm ||Z||^2 = sum(sum(A .* A .* W)).   
%
%   [D, X, Y, ITS] = SDDWEIGHT(...) also returns the number of inner
%   iterations for each outer iteration.   
%
%   [D, X, Y, ITS, RHO] = SDDWEIGHT(...) also returns the norm-squared of
%   the residual after each outer iteration. 
%
%   [D, X, Y, ITS, RHO, IITS] = SDDWEIGHT(...) also returns the number of
%   extra matrix-vector multiplies used in the initialization of each inner
%   iteration when using C=3 (see below).  
%
%   [...] = SDDWEIGHT(A, W, K) produces a rank-k SDDWEIGHT. The default
%   rank is 100. 
%
%   [...] = SDDWEIGHT(A, W, K, TOL) stops the inner iterations after the
%   improvement is less than TOL. The default value is 0.01.
%
%   [...] = SDDWEIGHT(A, W, K, TOL, L) specifies that the maximum number of
%   inner iterations is L. The default is 100.  
%
%   [...] = SDDWEIGHT(A, W, K, TOL, L, R) produces an SDDWEIGHT
%   approximation that is either of rank K or such that 
%   || A - X * diag(D) * Y' ||_W < R. The default is zero.
%
%   [...] = SDDWEIGHT(A, W, K, TOL, L, R, C) sets the choice for
%   initializing y in the inner iterations as follows:  
%
%      C = 1       Threshold
%      C = 2       Cycling
%      C = 3,      All elements of y are set to 1.
%      C = 4,      Every 100th element of y is set to 1 starting with 1.
%
%   Default is C = 1.
%
%
%SDDPACK: Software for the Semidiscrete Decomposition.
%Copyright (c) 1999 Tamara G. Kolda and Dianne P. O'Leary. 

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the Free 
% Software Foundation; either version 2 of the License, or (at your option)
% any later version.  
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
% or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
% for more details.  
%
% You should have received a copy of the GNU General Public License along
% with this program; if not, write to the Free Software Foundation, Inc., 59
% Temple Place - Suite 330, Boston, MA 02111-1307, USA.  

%%% Check Input Arguments

if ~exist('A') 
  error('Incorrect number of inputs.');
end

if ~exist('kmax')
  kmax = 100;
end

if ~exist('alphamin')
  alphamin = 0.01;
end

if ~exist('lmax')
  lmax = 100;
end

if ~exist('rhomin')
  rhomin = 0;
else
  rhomin = rhomin^2;
end

if ~exist('yinit')
  yinit = 1;
end
idx = 1;                                % only used for yinit = 2

%%% Initialization

[m, n] = size(A);                       % size of A
rho = sum(sum(A .* A .* W));            % squared residual norm 
iitssav = [];

%%% Outer Loop

for k = 1 : kmax
  
  %%% Initialize y for Inner Loop

  switch yinit(1),                      % All Ones
    case 1,                             % Threshold
      snorm = 0;
      iits = 0;
      while snorm < (rho / n)
        y = zeros(n, 1); 
        y(idx) = 1;
        snorm = sum((A(:,idx) .* A(:,idx) .* W(:,idx)));
        idx = mod(idx, n) + 1;
        iits = iits + 1;
      end
      iitssav(k) = iits;
    case 2,                             % Cycling Periodic Ones
      y = zeros(n, 1);
      y(mod(k-1,n)+1) = 1;
    case 3,
      y = ones(n, 1);
    case 4,                             % Periodic Ones
      y = zeros(n, 1);
      for i = 1 : 100 : n
        y(i) = 1;
      end
    otherwise,
      error('Invalid choice for C.');
  end % switch on yinit

  %%% Inner Loop

  for l = 1 : lmax

    %%% Fix y and Solve for x

    s = (A .* W) * y;
    v = W *(y .* y);

    [x, xcnt] = solveweight(s, v, m);
      
    %%% Fix x and Solve for y

    s = (A .* W)' * x;
    v = W' * (x .* x);
    
    [y, ycnt] = solveweight(s, v, n);
    
    %%% Check Progress
    
    beta = (s' * y)^2 / (v' * (y .* y));
    
    if (l > 1)
      alpha = (beta - betabar) / betabar;
      if (alpha <= alphamin)
        break
      end
    end

    betabar = beta;

  
  end % l-loop
  
  %%% Save
    
  d = (s' * y) / (v' * (y .* y));
  xsav(:, k) = x;
  ysav(:, k) = y;
  dsav(k, 1) = d;
  A = A - x * d * y';
  rho = rho - beta;
  rhosav(k) = rho;
  itssav(k) = l;
  
  %%% Threshold Test

  if (rho < rhomin)
    break;
  end
  
end % k-loop

return

%----------------------------------------------------------------------%

function [x, imax] = solveweight(s, v, m)

%SOLVEWEIGHT Solve SDDWEIGHT subproblem
%   [X] = SOLVE(S, V, M) computes max(X'S)/(X.*X)'*V where M is the size of
%   S.  
%
%   [X, I] = SOLVE(S, V, M) additionally returns number of nonzeros in X.
%
%For use with SDDWEIGHT.
%Tamara G. Kolda, Oak Ridge National Laboratory, 1999.
%Dianne P. O'Leary, University of Maryland and ETH, 1999.

if s == 0, error('s is zero'), end;

for i = 1 : m
  if s(i) < 0
    x(i, 1) = -1;
    s(i) = -s(i);
    ratio(i) = -s(i)/v(i);
  elseif s(i) > 0
    x(i, 1) = 1;
    ratio(i) = -s(i)/v(i);
  else
    x(i, 1) = 0;
    ratio(i) = 0;
  end
end

[sorts, indexsort] = sort(ratio);
clear f
num = s(indexsort(1));
den = v(indexsort(1));
f(1) = num^2/den;
for i = 2 : m
  num = num + s(indexsort(i));
  den = den + v(indexsort(i));
  f(i) = num^2/den;
end

imax = 1;
fmax = f(1);
for i = 2 : m
  if f(i) >= fmax
    imax = i;
    fmax = f(i);
  end
end

for i = (imax + 1) : m
  x(indexsort(i)) = 0;
end

return




