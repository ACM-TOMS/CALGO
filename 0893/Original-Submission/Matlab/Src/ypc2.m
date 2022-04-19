function [yp,ier] = ypc2(x,y,sigma,isl1,isln,bv1,bvn)
% ypc2:  C^2 global derivative estimates
%
% USAGE:  [yp,ier] = ypc2(x,y,sigma,isl1,isln,bv1,bvn);
%
%   This function solves a linear system for a set of
% first derivatives YP associated with a Hermite interpola-
% tory tension spline H(x).  The derivatives are chosen so
% that H(x) has two continuous derivatives for all x and H
% satisfies user-specified end conditions.
%
% On input:
%
%       X = Vector of length N containing a strictly in-
%           creasing sequence of abscissae:  X(I) < X(I+1)
%           for I = 1 to N-1.  N >= 2.
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
%       ISL1 = Option indicator for the condition at X(1):
%              ISL1 = 0 if YP(1) is to be estimated inter-
%                       nally by a constrained parabolic
%                       fit to the first three points.
%                       This is identical to the method used
%                       by Function YPC1.  BV1 is not used
%                       in this case.
%              ISL1 = 1 if the first derivative of H at X(1)
%                       is specified by BV1.
%              ISL1 = 2 if the second derivative of H at
%                       X(1) is specified by BV1.
%              ISL1 = 3 if YP(1) is to be estimated inter-
%                       nally from the derivative of the
%                       tension spline (using SIGMA(1))
%                       that interpolates the first three
%                       data points and has third derivative
%                       equal to zero at X(1).  Refer to
%                       ENDSLP.  BV1 is not used in this
%                       case.
%
%       ISLN = Option indicator for the condition at X(N):
%              ISLN = 0 if YP(N) is to be estimated inter-
%                       nally by a constrained parabolic
%                       fit to the last three data points.
%                       This is identical to the method used
%                       by Function YPC1.  BVN is not used
%                       in this case.
%              ISLN = 1 if the first derivative of H at X(N)
%                       is specified by BVN.
%              ISLN = 2 if the second derivative of H at
%                       X(N) is specified by BVN.
%              ISLN = 3 if YP(N) is to be estimated inter-
%                       nally from the derivative of the
%                       tension spline (using SIGMA(N-1))
%                       that interpolates the last three
%                       data points and has third derivative
%                       equal to zero at X(N).  Refer to
%                       ENDSLP.  BVN is not used in this
%                       case.
%
%       BV1,BVN = Boundary values or dummy parameters as
%                 defined by ISL1 and ISLN.
%
% On output:
%
%       YP = Vector of size(X) containing derivatives  of
%            H at the abscissae.  YP is zeros if IER ~= 0.
%
%       IER = Error indicator:
%             IER = 0 if no errors were encountered.
%             IER = 1 if N, ISL1, or ISLN is outside its
%                     valid range.
%             IER = I if X(I) <= X(I-1) for some I in the
%                     range 2 to N.
%
% Modules required by YPC2:  ENDSLP, SNHCSH, YPCOEF
%
%***********************************************************

m = size(x);
n = length(x);

yp = zeros(m);

if (n < 2  ||  isl1 < 0  ||  isl1 > 3  ||  isln < 0  || ...
    isln > 3)
   ier = 1;
   return
end
nm1 = n - 1;

% Set YP1 and YPN to the endpoint values.

if (isl1 == 0)
   if (n > 2)
      yp1 = endslp(x(1),x(2),x(3),y(1),y(2),y(3),0.0);
   end
elseif (isl1 ~= 3) 
   yp1 = bv1;
else
   if (n > 2) 
      yp1 = endslp(x(1),x(2),x(3),y(1),y(2),y(3),sigma(1));
   end
end
if (isln == 0) 
   if (n > 2) 
      ypn = endslp(x(n),x(nm1),x(n-2),y(n),y(nm1),y(n-2),0.0);
   end
elseif (isln ~= 3)
   ypn = bvn;
else
   if (n > 2) 
      ypn = endslp(x(n),x(nm1),x(n-2),y(n),y(nm1),y(n-2),sigma(nm1));
   end
end

% Arrays:

dx = diff(x);
if (any(dx <= 0))
   ier = find(dx <= 0, 1) + 1;
   return;
end
s = diff(y)./dx;
sig = abs(sigma);
wk = zeros(nm1,1);

% Solve the symmetric positive-definite tridiagonal linear
%   system.  The forward elimination step consists of div-
%   iding each row by its diagonal entry, then introducing a
%   zero below the diagonal.  This requires saving only the
%   superdiagonal (in WK) and the right hand side (in YP).

if (n == 2) 
   if (isl1 == 0  ||  isl1 == 3), yp1 = s(1); end
   if (isln == 0  ||  isln == 3), ypn = s(1); end
end

% Begin forward elimination.

[d1,sd1] = ypcoef(sig(1),dx(1));
r1 = (sd1+d1)*s(1);
wk(1) = 0;
yp(1) = yp1;
if (isl1 == 2) 
   wk(1) = sd1/d1;
   yp(1) = (r1-yp1)/d1;
end
for i = 2:nm1
   [d2,sd2] = ypcoef(sig(i),dx(i));
   r2 = (sd2+d2)*s(i);
   d = d1 + d2 - sd1*wk(i-1);
   wk(i) = sd2/d;
   yp(i) = (r1 + r2 - sd1*yp(i-1))/d;
   d1 = d2;
   sd1 = sd2;
   r1 = r2;
end
d = d1 - sd1*wk(nm1);
yp(n) = ypn;
if (isln == 2) 
   yp(n) = (r1 + ypn - sd1*yp(nm1))/d;
end

% Back substitution:

for i = nm1:-1:1
   yp(i) = yp(i) - wk(i)*yp(i+1);
end
ier = 0;
return;

end  % ypc2
