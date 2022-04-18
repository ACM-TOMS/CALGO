function [yp,ier] = ypc1(x,y)
% ypc1:  Local derivative estimates
%
% USAGE:  [yp,ier] = ypc1(x,y);
%
%   This function employs a three-point quadratic interpo-
% lation method to compute local derivative estimates YP
% associated with a set of N data points.  The interpolation
% formula is the monotonicity-constrained parabolic method
% described in the reference cited below.  A Hermite int-
% erpolant of the data values and derivative estimates pre-
% serves monotonicity of the data.  Linear interpolation is
% used if N = 2.  The method is invariant under a linear
% scaling of the coordinates but is not additive.
%
% On input:
%
%       X = Vector of length N containing a strictly in-
%           creasing sequence of abscissae:  X(I) < X(I+1)
%           for I = 1 to N-1.  N >= 2.
%
%       Y = Vector of length N containing data values asso-
%           ciated with the abscissae.
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
%             IER = 1 if N < 2.
%             IER = I if X(I) <= X(I-1) for some I in the
%                     range 2 to N.
%
% Reference:  J. M. Hyman, "Accurate Monotonicity-preserving
%               Cubic Interpolation",  LA-8796-MS, Los
%               Alamos National Lab, Feb. 1982.
%
% Modules required by YPC1:  None
%
%***********************************************************

n = length(x);

yp = zeros(size(x));

if (n < 2)
   ier = 1;
   return;
end

dx = diff(x);
if (any(dx <= 0))
   ier = find(dx <= 0, 1) + 1;
   return;
end
s = diff(y)./dx;
if (n == 2) 

% Use linear interpolation for N = 2.

   yp(1) = s;
   yp(2) = s;
   ier = 0;
   return;
end

% N >= 3.  YP(1) = S(1) + DX(1)*(S(1)-S(2))/(DX(1)+DX(2)) 
%   unless this results in YP(1)*S(1) <= 0 or abs(YP(1)) > 
%   3*abs(S(1)).

t = s(1) + dx(1)*(s(1)-s(2))/(dx(1)+dx(2));
if (s(1) >= 0) 
   yp(1) = min([max([0, t]), 3.0*s(1)]);
else
   yp(1) = max([min([0, t]), 3.0*s(1)]);
end

% For i = 2 to N-1, YP(i) = (DX(i-1)*S(i)+DX(i)*S(i-1))/
%   (DX(i-1)+DX(i)) subject to the constraint that YP(i) 
%   has the sign of either S(i-1) or S(i), whichever has 
%   larger magnitude, and abs(YP(i)) <= 3*min(abs(S(i-1)),
%   abs(S(i))).

smag = abs(s);
i = 2:n-1;
yp(i) = (dx(i-1).*s(i) + dx(i).*s(i-1))./(dx(i-1)+dx(i));

% Partition the interior point indices i = 2:n-1 into subsets
%   k1 and k2.

k1 = find(smag(i-1) > smag(i)) + 1;
k2 = setdiff(i, k1);

yp(k1(sign(yp(k1)).*sign(s(k1-1)) < 0)) = 0;
k = find(abs(yp(k1)) > 3.0*smag(k1));
yp(k1(k)) = 3.0*smag(k1(k)).*sign(yp(k1(k)));

yp(k2(sign(yp(k2)).*sign(s(k2)) <= 0)) = 0;
k = find(abs(yp(k2)) > 3.0*smag(k2-1));
yp(k2(k)) = 3.0*smag(k2(k)-1).*sign(yp(k2(k)));

% YP(N) = S(N-1) + DX(N-1)*(S(N-1)-S(N-2))/(DX(N-2)+DX(N-1))
%   subject to the constraint that YP(N) has the sign of 
%   S(N-1) and abs(YP(N)) <= 3*abs(S(N-1)).

t = s(n-1) + dx(n-1)*(s(n-1)-s(n-2))/(dx(n-2)+dx(n-1));
if (s(n-1) >= 0) 
   yp(n) = min([max([0, t]), 3.0*s(n-1)]);
else
   yp(n) = max([min([0, t]), 3.0*s(n-1)]);
end

% No error encountered.

ier = 0;
return;

end  % ypc1
