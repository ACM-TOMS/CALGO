%
%   File:         vcalls2.m
%
%   Description:  Example of how to define struct parameters
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

name = '(I-2uu'') D (I-2uu'')';

uutmatvecpar.d = rand(50,1);
uutmatvecpar.u = rand(50,1);
uutmatvecpar.u = uutmatvecpar.u/norm(uutmatvecpar.u);

g     = rand(50,1);
Delta = 1;

%
% These statements redefine the values of epsilon.Delta, epsilon.HC
% lstrs will use the new values
%

epsilon.Delta = 1e-3;
epsilon.HC    = 1e-8;

[x,lambda,info,moreinfo] = lstrs(@uutmatvec,g,Delta,epsilon,[],[],uutmatvecpar);

%
% These statements redefine the values of opts.tol, opts.p opts.v0 for eigs_lstrs
% and, lopts.message_level and lopts.name
%
epar.tol = 1e-3;
epar.p   = 15;
epar.v0  = ones(51,1)/sqrt(51); % Initial vector for ARPACK

lopts.message_level = 2;
lopts.name          = name;

%[x,lambda,info,moreinfo] =  ...
%lstrs(@uutmatvec,g,Delta,epsilon,[],lopts,uutmatvecpar,epar);

% The longest possible call to lstrs

%[x,lambda,info,moreinfo] =  ...
%lstrs(@uutmatvec,g,Delta,epsilon,@tcheigs_lstrs_gateway,lopts,uutmatvecpar,epar);

