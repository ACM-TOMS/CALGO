function [sinhm,coshm,coshmm] = snhcsh(x)
% snhcsh:  Modified hyperbolic function evaluation
%
% USAGE:  [sinhm,coshm,coshmm] = snhcsh(x);
%
%   This function evaluates the modified hyperbolic 
% functions sinh(x)-x, cosh(x)-1, and cosh(x)-1-x/2 with
% machine precision accuracy (relative error bounded by 
% 3.4E-20 for a floating point number system with
% sufficient precision).
%
% On input:
%
%       X = Point or vector of points at which the functions
%           are to be evaluated.
%
% On output:
%
%       SINHM = sinh(X) - X.
%
%       COSHM = cosh(X) - 1.
%
%       COSHMM = cosh(X) - 1 - X*X/2.
%
% Modules required by SNHCSH:  None
%
%***********************************************************

% Coefficients defining rational approximations for small x.

p1 = -3.51754964808151394800e5;
p2 = -1.15614435765005216044e4;
p3 = -1.63725857525983828727e2;
p4 = -7.89474443963537015605e-1;
q1 = -2.11052978884890840399e6;
q2 = 3.61578279834431989373e4;
q3 = -2.77711081420602794433e2;
q4 = 1.0;

ax = abs(x);
xs = ax.*ax;

m = size(x);
sinhm = zeros(m);
coshm = zeros(m);
coshmm = zeros(m);

% Approximations for small X:

k = find(ax <= 0.5);
xc = x(k).*xs(k);
p = ((p4*xs(k)+p3).*xs(k)+p2).*xs(k)+p1;
q = ((q4*xs(k)+q3).*xs(k)+q2).*xs(k)+q1;
sinhm(k) = xc.*(p./q);
xsd4 = 0.25*xs(k);
xsd2 = xsd4 + xsd4;
p = ((p4*xsd4+p3).*xsd4+p2).*xsd4+p1;
q = ((q4*xsd4+q3).*xsd4+q2).*xsd4+q1;
f = xsd4.*(p./q);
coshmm(k) = xsd2.*f.*(f+2.0);
coshm(k) = coshmm(k) + xsd2;

% Approximations for large X:

k = find(ax > 0.5);
expx = exp(ax(k));
sinhm(k) = -(((1.0./expx+ax(k))+ax(k))-expx)/2.0;
k1 = find(x < 0);
k1 = intersect(k,k1);
sinhm(k1) = -sinhm(k1);
coshm(k) = ((1.0./expx-2.0)+expx)/2.0;
coshmm(k) = coshm(k) - xs(k)/2.0;
return;

end  % snhcsh
