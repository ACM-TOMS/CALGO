%
%   File:         vcalls1.m
%
%   Description:  Example of valid calls to lstrs
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

name   ='Identity';
H      = eye(50);
g      = ones(50,1);
mu     = -3;
xexact = -ones(50,1)/(1-mu);
Delta  = norm(xexact);

% Eigensolver is tcheigs_lstrs_gateway
% Initial vector for ARPACK is random

[x,lambda,info,moreinfo] = lstrs(@mv,g,Delta,[],@tcheigs_lstrs_gateway);
%[x,lambda,info,moreinfo] = lstrs(@mv,g,Delta,[],'tcheigs_lstrs_gateway');

% Eigensolver is eig_gateway

%[x,lambda,info,moreinfo] = lstrs(H,g,Delta,[],@eig_gateway);

%
% Defining maxiter, message_level, name
% Default values are used for the remaining parameters
%
lopts.maxiter       = 3;
lopts.message_level = 0;
lopts.name          = name;

%[x,lambda,info,moreinfo] = lstrs(@mv,g,Delta,[],[],lopts);

