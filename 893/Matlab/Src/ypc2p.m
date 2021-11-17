function [yp,ier] = ypc2p(x,y,sigma)
% ypc2p:  C^2 global derivative estimates, periodic case
%
% USAGE:  [yp,ier] = ypc2p(x,y,sigma);
%
%   This function solves a linear system for a set of
% first derivatives YP associated with a Hermite interpola-
% tory tension spline H(x).  The derivatives are chosen so
% that H(x) has two continuous derivatives for all x, and H
% satisfies periodic end conditions:  first and second der-
% ivatives of H at X(1) agree with those at X(N), and thus
% the length of a period is X(N) - X(1).  It is assumed that
% Y(N) = Y(1).
%
% On input:
%
%       X = Vector of length N containing a strictly in-
%           creasing sequence of abscissae:  X(I) < X(I+1)
%           for I = 1 to N-1.  N >= 3.
%
%       Y = Vector of length N containing data values asso-
%           ciated with the abscissae.  H(X(I)) = Y(I) for
%           I = 1 to N.
%
%       SIGMA = Vector of length N-1 containing tension
%               factors.  SIGMA(I) is associated with inter-
%               val (X(I),X(I+1)) for I = 1 to N-1.  If
%               SIGMA(I) = 0, H is the Hermite cubic interp-
%               olant of the data values and computed deriv-
%               atives at X(I) and X(I+1), and if all
%               tension factors are zero, H is the C-2 cubic
%               spline interpolant that satisfies the end
%               conditions.
%
% On output:
%
%       YP = Vector of size(X) containing derivatives of
%            H at the abscissae.  YP is zeros if IER ~= 0.
%
%       IER = Error indicator:
%             IER = 0 if no errors were encountered.
%             IER = 1 if N is outside its valid range.
%             IER = I if X(I) <= X(I-1) for some I in the
%                     range 2 to N.
%
% Modules required by YPC2P:  SNHCSH, YPCOEF
%
%***********************************************************

m = size(x);
n = length(x);

yp = zeros(m);

if (n < 3)
   ier = 1;
   return;
end
nm1 = n - 1;
nm2 = n - 2;
nm3 = n - 3;
np1 = n + 1;

% Arrays:

dx = diff(x);
if (any(dx <= 0))
   ier = find(dx <= 0, 1) + 1;
   return;
end
s = diff(y)./dx;
sig = abs(sigma);
wk = zeros(2*nm1,1);

% The system is order N-1, symmetric, positive-definite, and
%   tridiagonal except for nonzero elements in the upper
%   right and lower left corners.  The forward elimination
%   step zeros the subdiagonal and divides each row by its
%   diagonal entry for the first N-2 rows.  The superdiago-
%   nal is stored in WK(I), the negative of the last column
%   (fill-in) in WK(N+I), and the right hand side in YP(I)
%   for I = 1 to N-2.

% i = nm1

[dnm1,sdnm1] = ypcoef(sig(nm1),dx(nm1));
rnm1 = (sdnm1+dnm1)*s(nm1);

% i = 1

[d1,sd1] = ypcoef (sig(1),dx(1));
r1 = (sd1+d1)*s(1);
d = dnm1 + d1;
wk(1) = sd1/d;
wk(np1) = -sdnm1/d;
yp(1) = (rnm1+r1)/d;
for i = 2:nm2
   [d2,sd2] = ypcoef(sig(i),dx(i));
   r2 = (sd2+d2)*s(i);
   d = d1 + d2 - sd1*wk(i-1);
   din = 1.0/d;
   wk(i) = sd2*din;
   npi = n + i;
   wk(npi) = -sd1*wk(npi-1)*din;
   yp(i) = (r1 + r2 - sd1*yp(i-1))*din;
   sd1 = sd2;
   d1 = d2;
   r1 = r2;
end

% The backward elimination step zeros the superdiagonal
%   (first N-3 elements).  WK(I) and YP(I) are overwritten
%   with the negative of the last column and the new right
%   hand side, respectively, for I = N-2, N-3, ..., 1.

npi = n + nm2;
wk(nm2) = wk(npi) - wk(nm2);
for i = nm3:-1:1
   yp(i) = yp(i) - wk(i)*yp(i+1);
   npi = n + i;
   wk(i) = wk(npi) - wk(i)*wk(i+1);
end

% Solve the last equation for YP(N-1).

ypnm1 = (r1 + rnm1 - sdnm1*yp(1) - sd1*yp(nm2))/ ...
        (d1 + dnm1 + sdnm1*wk(1) + sd1*wk(nm2));

% Back substitute for the remainder of the solution
%   components.

yp(nm1) = ypnm1;
for i = 1:nm2
   yp(i) = yp(i) + wk(i)*ypnm1;
end

% YP(N) = YP(1).

yp(n) = yp(1);
ier = 0;
return;

end  % ypc2p
