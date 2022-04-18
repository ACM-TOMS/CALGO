function [nqpts,maxturn,maxlen,maxrefn] = scplotopt(options)
%SCPLOTOPT Parameters used by S-C plotting routines.
%       OPTIONS(1): Number of quadrature points per integration.
%                   Approximately equals -log10(error).  Increase if plot
%                   has false little zigzags in curves (default 4). 
%       OPTIONS(2): Maximum allowed turning angle at each plotted point,
%                   in degrees (default 12).
%       OPTIONS(3): Max allowed line segment length, as a proportion of the
%                   largest finite polygon side (default 0.05).
%       OPTIONS(4): Max allowed number of adaptive refinements made to meet
%                   other requirements (default 10).
%      
%       See also HPPLOT, DPLOT, DEPLOT, STPLOT, RPLOT.
%
%	Written by Toby Driscoll.  Last updated 5/24/95.

user = options;
lenu = length(user);
options = zeros(1,4);
options(1:lenu) = user(1:lenu);
options = options + (options==0).*[4,12,.05,10];

nqpts = options(1);
maxturn = options(2);
maxlen = options(3);
maxrefn = options(4);

