function [dsav, xsav, ysav, itssav, rhosav, iitssav] = ...
    sdd(A, kmax, alphamin, lmax, rhomin, yinit)

%SDD  Semidiscrete Decomposition.
%
%   [D, X, Y] = SDD(A) produces discrete matrices X and Y and a vector D
%   such that X * diag(D) * Y' is the 100-term SDD that approximates A.
%
%   [D, X, Y, ITS] = SDD(...) also returns the number of inner iterations
%   for each outer iteration.  
%
%   [D, X, Y, ITS, RHO] = SDD(...) also returns a vector RHO containing the
%   norm-squared of the residual after each outer iteration.
%
%   [D, X, Y, ITS, RHO, IITS] = SDD(...) also returns a vector IITS
%   containing the number of extra matrix-vector multiplies used in the
%   initialization of each inner iteration when using C=1 (see below). 
%
%   [...] = SDD(A, K) produces a K-term SDD. The default is 100.
%
%   [...] = SDD(A, K, TOL) stops the inner iterations after the
%   improvement is less than TOL. The default is 0.01.
%
%   [...] = SDD(A, K, TOL, L) specifies that the maximum number of inner
%   iterations is L. The default is 100.
%
%   [...] = SDD(A, K, TOL, L, R) produces an SDD approximation that
%   either has K terms or such that norm(A - X * diag(D) * Y', 'fro') < R.
%   The default is zero.
%
%   [...] = SDD(A, K, TOL, L, R, C) sets the choice for initializing y in
%   the inner iterations as follows: 
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
  error('Incorrent number of inputs.');
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
idx = 1;				% only used for yinit = 1

%%% Initialization

[m, n] = size(A);			% size of A
rho = norm(A, 'fro')^2;			% squared residual norm 
iitssav = [];

%%% Outer Loop

for k = 1 : kmax
  
  %%% Initialize y for Inner Loop

  switch yinit,				
    case 1,				% Threshold
      s = zeros(m, 1);
      iits = 0;
      while (norm(s)^2) < (rho / n)
	y = zeros(n, 1); 
	y(idx) = 1;
	s = A * y;
	if k > 1
	  s = s - (xsav * (dsav .* (ysav' * y)));
	end
	idx = mod(idx, n) + 1;
	iits = iits + 1;
      end
      iitssav(k) = iits;
    case 2,				% Cycling Periodic Ones
      y = zeros(n, 1);
      y(mod(k-1,n)+1) = 1;
    case 3,				% All Ones
      y = ones(n, 1);
    case 4,				% Periodic Ones
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

    s = A * y;
    if k > 1
      s = s - (xsav * (dsav .* (ysav' * y)));
    end

    [x, xcnt] = sddsolve(s, m);
      
    %%% Fix x and Solve for y

    s = A' * x;
    if (k > 1)
      s = s - (ysav * (dsav .* (xsav' * x)));
    end
    
    [y, ycnt, fmax] = sddsolve(s, n);
    
    %%% Check Progress
    
    d = sqrt(fmax * ycnt) / (ycnt * xcnt);
      
    beta = d^2 * ycnt * xcnt;
    
    if (l > 1)
      alpha = (beta - betabar) / betabar;
      if (alpha <= alphamin)
	break
      end
    end

    betabar = beta;
  
  end % l-loop
  
  %%% Save
    
  xsav(:, k) = x;
  ysav(:, k) = y;
  dsav(k, 1) = d;
  rho = max([rho - beta, 0]);
  rhosav(k) = rho;
  itssav(k) = l;
  
  %%% Threshold Test

  if (rho <= rhomin)
    break;
  end
  
end % k-loop

return

%----------------------------------------------------------------------%

function [x, imax, fmax] = sddsolve(s, m)

%SDDSOLVE Solve SDD subproblem
%
%   [X] = SDDSOLVE(S, M) computes max (X' * S) / (X' * X) where M is the
%   size of S.  
%
%   [X, I] = SDDSOLVE(S, M) additionally returns number of nonzeros in X.
%
%   [X, I, F] = SDDSOLVE(S, M) additionally returns value of function at the
%   optimum.  
%
%For use with SDD.
%Tamara G. Kolda, Oak Ridge National Laboratory, 1999.
%Dianne P. O'Leary, University of Maryland and ETH, 1999.

for i = 1 : m
  if s(i) < 0
    x(i, 1) = -1;
    s(i) = -s(i);
  else
    x(i, 1) = 1;
  end
end

[sorts, indexsort] = sort(-s);
sorts = -sorts;

clear f
f(1) = sorts(1);
for i = 2 : m
  f(i) = sorts(i) + f(i - 1);
end

f = (f.^2) ./ [1:m];

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


