function [yp,ier] = ypc1p(x,y)
% ypc1p:  Local derivative estimates, periodic case
%
% USAGE:  [yp,ier] = ypc1p(x,y);
%
%   This function employs a three-point quadratic interpo-
% lation method to compute local derivative estimates YP
% associated with a set of N data points (X(I),Y(I)).  It
% is assumed that Y(N) = Y(1), and YP(N) = YP(1) on output.
% Thus, a Hermite interpolant H(x) defined by the data
% points and derivative estimates is periodic with period
% X(N)-X(1).  The derivative-estimation formula is the
% monotonicity-constrained parabolic fit described in the
% reference cited below:  H(x) is monotonic in intervals in
% which the data is monotonic.  The method is invariant
% under a linear scaling of the coordinates but is not
% additive.
%
% On input:
%
%       X = Vector of length N containing a strictly in-
%           creasing sequence of abscissae:  X(I) < X(I+1)
%           for I = 1 to N-1.  N >= 3.
%
%       Y = Vector of length N containing data values asso-
%           ciated with the abscissae.  Y(N) = Y(1).
%
% On output:
%
%       YP = Vector of size(X) containing estimated deriv-
%            atives at the abscissae unless IER ~= 0.
%            YP is zeros if IER = 1, and is only partially
%            computed if IER > 1.
%
%       IER = Error indicator:
%             IER = 0 if no errors were encountered.
%             IER = 1 if N < 3 or Y(N) ~= Y(1).
%             IER = I if X(I) <= X(I-1) for some I in the
%                     range 2 to N.
%
% Reference:  J. M. Hyman, "Accurate Monotonicity-preserving
%               Cubic Interpolation",  LA-8796-MS, Los
%               Alamos National Lab, Feb. 1982.
%
% Modules required by YPC1P:  None
%
%***********************************************************

n = length(x);

yp = zeros(size(x));

if (n < 3  ||  y(n) ~= y(1))
   ier = 1;
   return;
end

% For i = 2 to N-1, YP(i) = (DX(i-1)*S(i)+DX(i)*S(i-1))/
%   (DX(i-1)+DX(i)) subject to the constraint that YP(i) 
%   has the sign of either S(i-1) or S(i), whichever has 
%   larger magnitude, and abs(YP(i)) <= 3*min(abs(S(i-1)),
%   abs(S(i))).

dx = diff(x);
if (any(dx <= 0))
   ier = find(dx <= 0, 1) + 1;
   return;
end
s = diff(y)./dx;
smag = abs(s);
i = 2:n-1;
yp(i) = (dx(i-1).*s(i) + dx(i).*s(i-1))./(dx(i-1)+dx(i));

% Partition the interior point indices 2:n-1 into subsets
%   k1 and k2.

k1 = find(smag(i-1) > smag(i)) + 1;
k2 = setdiff(i, k1);

yp(k1(sign(yp(k1)).*sign(s(k1-1)) < 0)) = 0;
k = find(abs(yp(k1)) > 3.0*smag(k1));
yp(k1(k)) = 3.0*smag(k1(k)).*sign(yp(k1(k)));

yp(k2(sign(yp(k2)).*sign(s(k2)) <= 0)) = 0;
k = find(abs(yp(k2)) > 3.0*smag(k2-1));
yp(k2(k)) = 3.0*smag(k2(k)-1).*sign(yp(k2(k)));

% Compute YP(N) = YP(1):  i = 1 and i-1 = n-1.

t = (dx(n-1)*s(1) + dx(1)*s(n-1))/(dx(n-1)+dx(1));
if (smag(n-1) > smag(1))
   sgn = sign(s(n-1));
else
   sgn = sign(s(1));
end
if (sgn > 0) 
   yp(1) = min([max([0, t]),  3.0*min([smag(n-1), smag(1)])]);
else
   yp(1) = max([min([0, t]), -3.0*min([smag(n-1), smag(1)])]);
end
yp(n) = yp(1);

% No error encountered.

ier = 0;
return;

end  % ypc1p
