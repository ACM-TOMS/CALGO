%
%   File:         regularization.m
%
%   Description:  Solving a regularization problem with lstrs.
%                 Problem 'phillips' from Regularization Tools 
%                 by Per Christian Hansen
%
%   Authors:      Marielba Rojas
%                    mr@imm.dtu.dk
%                 Sandra Santos
%                 Danny Sorensen
%
%   Version:      1.2
%
%   System:       MATLAB 6.0 or higher 
%
%   Date:         15 March 2007
%

%
% Note: 
% problem 'phillips' is part of the Regularization Tools Package
% by Per Christian Hansen, available from
% http://www.imm.dtu.dk/~pch/Regutools/index.html
%

n = 300;
[A,b,xexact] = phillips(n);

atamvpar.A = A;
g          = - (b'*A)';
Delta      = norm(xexact); 

lopts.name       = 'phillips';
lopts.plot       = 'y';
lopts.correction = 'n';
lopts.interior   = 'n';

epsilon.Delta = 1e-2;

% initial vector for first call to eigs in lstrs
epar.v0 = ones(n+1,1)/sqrt(n+1);

%
% This yields a quasioptimal solution
%
[x,lambda,info,moreinfo] = ...
lstrs(@atamv,g,Delta,epsilon,@tcheigs_lstrs_gateway,lopts,atamvpar,epar);


% The following yields a Boundary Solution that is closer to
% the exact solution and uses a lower number of iterations.
% The eigensolver is the default: eigs_lstrs_gateway.
% Uncomment to try

%lopts.heuristics = 1;
%[x,lambda,info,moreinfo] = ...
%lstrs(@atamv,g,Delta,epsilon,[],lopts,atamvpar,epar);

