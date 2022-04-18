function [ode,newton,tol,maxiter] = scimapopt(options)
%SCIMAPOPT Parameters used by S-C inverse-mapping routines.
%       OPTIONS(1): Algorithm (default 0)
%	            0--use ode to get initial guess, then Newton iters.
%	            1--use ode only
%	            2--use Newton only; take Z0 as initial guess
%       OPTIONS(2): Error tolerance for solution (default 1e-8)
%	OPTIONS(3): Maximum number of Newton iterations (default 10)
%
%	See also HPINVMAP, DINVMAP, DEINVMAP, RINVMAP, STINVMAP.
%	
%       Written by Toby Driscoll.  Last updated 5/24/95.

user = options;
lenu = length(user);
options = zeros(1,3);
options(1:lenu) = user(1:lenu);
options = options + (options==0).*[0,1e-8,10];

ode = options(1)==0 | options(1)==1;
newton = options(1)==0 | options(1)==2;
tol = options(2);
maxiter = options(3);

