function [y,d] = dfixwc(w,beta,z,c,wc,z0)
%DFIXWC Fix conformal center of disk map.
%       The conformal center WC of a Schwarz-Christoffel interior disk
%       map is defined as the image of zero.  The parameter problem
%       solver DPARAM does not allow control over the placement of the
%       conformal center.  Using the output Z,C from DPARAM, [Z0,C0] =
%       DFIXWC(W,BETA,Z,C,WC) computes a Moebius transformation so that
%       if Z0 and C0 are used in place of Z and C, the conformal center
%       of the resulting map will be WC.
%
%	See also DPARAM, PTSOURCE.
%	
%	Written by Toby Driscoll.  Last updated 5/24/95.

n = length(w);

if nargin < 6
  z0 = [];
end 

zc = dinvmap(wc,w,beta,z,c,[],z0,[0,1e-10]);

% Transform prevertices.
y = ((1-zc')/(1-zc))*(z-zc)./(1-zc'*z);
y(n) = 1;				% force it to be exact

% Recalculate constant from scratch.
mid = (y(1)+y(2))/2;
qdat = scqdata(beta,10);
d = (w(1) - w(2))/...
    (dquad(y(2),mid,2,y,beta,qdat) - dquad(y(1),mid,1,y,beta,qdat));

