%
%   File:         simple.m
%
%   Description:  A simple trust-region problem where the Hessian
%                 is the Identity matrix
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

name   = 'Identity';
H      = eye(50);
g      = ones(50,1);
mu     = -3;
xexact = -ones(50,1)/(1-mu);
Delta  = norm(xexact);

%
% The simplest possible calls to lstrs. Default values are used.
%

[x,lambda,info,moreinfo] = lstrs(H,g,Delta);
%[x,lambda,info,moreinfo] = lstrs(@mv,g,Delta);

