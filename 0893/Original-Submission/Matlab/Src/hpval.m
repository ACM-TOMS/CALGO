function [hp,ier] = hpval(t,x,y,yp,sigma)
% hpval:  First derivative of Hermite interpolatory tension spline
%
% USAGE:  [hp,ier] = hpval(t,x,y,yp,sigma);
%
%   This function evaluates the first derivative HP of a
% Hermite interpolatory tension spline H at one or more
% points
%
% On input:
%
%       T = Point or vector of points at which HP is to be 
%           evaluated.  Extrapolation is performed if T < 
%           X(1) or T > X(N).
%
%       X = Vector of length N containing the abscissae.
%           These must be in strictly increasing order:
%           X(I) < X(I+1) for I = 1 to N-1.  N >= 2.
%
%       Y = Vector of length N containing data values.
%           H(X(I)) = Y(I) for I = 1 to N.
%
%       YP = Vector of length N containing first deriva-
%            tives.  HP(X(I)) = YP(I) for I = 1 to N.
%
%       SIGMA = Vector of length N-1 containing tension fac-
%               tors whose absolute values determine the
%               balance between cubic and linear in each
%               interval.  SIGMA(I) is associated with int-
%               erval (I,I+1) for I = 1 to N-1.
%
% On output:
%
%       HP = Column vector of length(T) containing deriva-
%            tive values HP(T), or zeros if IER < 0.
%
%       IER = Optional error indicator:
%             IER = 0  if no errors were encountered and
%                      X(1) <= T <= X(N) for all components
%                      of T.
%             IER = 1  if no errors were encountered and
%                      extrapolation was necessary.
%             IER = -1 if the abscissae are not in strictly
%                      increasing order.  (This error will
%                      not necessarily be detected.)
%
% Module required by HPVAL:  SNHCSH
%
%***********************************************************

global SBIG

% Convert all input row vectors to column vectors.

t = t(:);
x = x(:);
y = y(:);
yp = yp(:);
sigma = sigma(:);

n = length(x);
m = size(t);
hp = zeros(m);

% Find the index vector I of the left endpoints of the intervals
%   containing the elements of T.  If T < X(1) or T > X(N), 
%   extrapolation is performed using the leftmost or rightmost 
%   interval.

i = ones(m);
for j = 2:n-1
   i(x(j) <= t) = j;
end
if (nargout > 1)
   ier = any(t < x(1)) || any(t > x(n));
end

% Compute interval widths DX, local coordinates B1 and B2,
%   and second differences D1 and D2.

ip1 = i + 1;
dx = x(ip1) - x(i);
if (nargout > 1  &&  any(dx <= 0))
   ier = -1;
   return;
end
b1 = (x(ip1) - t)./dx;
b2 = 1.0 - b1;
s1 = yp(i);
s = (y(ip1)-y(i))./dx;
d1 = s - s1;
d2 = yp(ip1) - s;
sig = abs(sigma(i));

% For SIG = 0, H is the Hermite cubic interpolant.

k = find(sig < 1.e-9);
hp(k) = s1(k) + b2(k).*(d1(k) + d2(k) - 3.0*b1(k).*(d2(k)-d1(k)));

% For 0 < SIG <= .5, use approximations designed to avoid
%   cancellation error in the hyperbolic functions.

k = find(sig >= 1.e-9  &  sig <= 0.5);
sb2 = sig(k).*b2(k);
[sm,cm,cmm] = snhcsh(sig(k));
[sm2,cm2] = snhcsh(sb2);
sinh2 = sm2 + sb2;
e = sig(k).*sm - cmm - cmm;
hp(k) = s1(k) + ((cm.*cm2-sm.*sinh2).*(d1(k)+d2(k)) + ...
                 sig(k).*(cm.*sinh2-(sm+sig(k)).*cm2).*d1(k))./e;

% For SIG > .5, use negative exponentials in order to avoid
%   overflow.  Note that EMS = EXP(-SIG).  In the case of
%   extrapolation (negative B1 or B2), H is approximated by
%   a linear function if -SIG*B1 or -SIG*B2 is large.

k = find(sig > 0.5);
sb1 = sig(k).*b1(k);
sb2 = sig(k) - sb1;
k1 = find(-sb1 > SBIG  |  -sb2 > SBIG);
k2 = k(k1);
hp(k2) = s(k2);
k1 = setdiff((1:length(k)),k1);
k = k(k1);
e1 = exp(-sb1(k1));
e2 = exp(-sb2(k1));
ems = e1.*e2;
tm = 1.0 - ems;
e = tm.*(sig(k).*(1.0+ems) - tm - tm);
hp(k) = s(k) + (tm.*((e2-e1).*(d1(k)+d2(k)) + tm.*(d1(k)-d2(k))) + ...
                sig(k).*((e1.*ems-e2).*d1(k) + (e1-e2.*ems).*d2(k)))./e;
return;

end  % hpval
