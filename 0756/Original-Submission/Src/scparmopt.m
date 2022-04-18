function [trace,tol] = scparmopt(options)
%SCPARMOPT Parameters used by S-C parameter problem routines.
%       OPTIONS(1): Nonzero causes some intermediate results to be
%                   displayed (default 0)
%       OPTIONS(2): Error tolerance for solution (default 1e-8)
%	
%	See also HPPARAM, DPARAM, DEPARAM, STPARAM, RPARAM.
%
%	Written by Toby Driscoll.  Last updated 5/24/95.

user = options;
lenu = length(user);
options = zeros(1,2);
options(1:lenu) = user(1:lenu);
options = options + (options==0).*[0,1e-8];

trace = options(1);
tol = options(2);

