function [ys,yp,ier] = smcrv(x,y,sigma,period,w,sm,smtol)
% smcrv:  C^2 smoothing curve
%
% USAGE:  [ys,yp,ier] = smcrv(x,y,sigma,period,w,sm,smtol);
%
%   Given a sequence of abscissae X with associated data
% values Y and tension factors SIGMA, this function deter-
% mines a set of function values YS and first derivatives YP
% associated with a Hermite interpolatory tension spline
% H that smoothes the data.  H has two continuous deriva-
% tives at every point, and satisfies either natural or
% periodic end conditions.  The values and derivatives are 
% chosen to minimize a quadratic functional Q1(YS,YP) 
% subject to the constraint Q2(YS) <= SM for Q2(YS) =
% (Y-YS)^T*W*(Y-YS), where ^T denotes transpose, and W is
% a diagonal matrix of positive weights.
%
%   Functions HVAL, HPVAL, HPPVAL, HPPVAL, and TSINTL may be
% called to compute values, derivatives, and integrals of H.  
% The function values YS must be used as data values in those
% functions.
%
%   The smoothing procedure is an extension of the method
% for cubic spline smoothing due to C. Reinsch:  Numer.
% Math., 10 (1967) and 16 (1971).  Q1 is defined as the sum
% of integrals over the intervals (X(I),X(I+1)) of HPP^2 +
% (SIGMA(I)/DX)^2*(HP-S)^2, where DX = X(I+1)-X(I), HP and
% HPP denote first and second derivatives of H, and S =
% (YS(I+1)-YS(I))/DX.  Introducing a smoothing parameter P,
% and assuming the constraint is active, the problem is
% equivalent to minimizing Q(P,YS,YP) = Q1(YS,YP) +
% P*(Q2(YS)-SM).  The procedure consists of finding a zero
% of G(P) = 1/SQRT(Q2) - 1/SQRT(SM), where YS and YP satisfy
% the order 2N symmetric positive-definite linear system
% obtained by setting the gradient of Q (treated as a func-
% tion of YS and YP) to zero.
%
%   Note that the interpolation problem corresponding to
% YS = Y, SM = 0, and P infinite is solved by Function
% YPC2 or YPC2P.
%
% On input:
%
%       X = Vector of length N containing a strictly in-
%           creasing sequence of abscissae:  X(I) < X(I+1)
%           for I = 1 to N-1.  N >= 2 if PERIOD = FALSE,
%           and N >= 3 if PERIOD = TRUE.
%
%       Y = Vector of length N containing data values assoc-
%           iated with the abscissae.  If PERIOD = TRUE, it
%           is assumed that Y(N) = Y(1).
%
%       SIGMA = Vector of length N-1 containing tension
%               factors.  SIGMA(I) is associated with inter-
%               val (X(I),X(I+1)) for I = 1 to N-1.  If
%               SIGMA(I) = 0, H is cubic, and as SIGMA in-
%               creases, H approaches linear in the inter-
%               val.
%
%       PERIOD = Periodic end condition flag:
%                PERIOD = 0 if H is to satisfy natural end
%                           conditions:  zero second der-
%                           ivatives at X(1) and X(N).
%                PERIOD = 1 if H is to satisfy periodic
%                           end conditions:  the values
%                           and first two derivatives at
%                           X(1) agree with those at X(N),
%                           and a period thus has length
%                           X(N)-X(1).
%
%       W = Vector of length N containing positive weights
%           associated with the data values.  The recommend-
%           ed value of W(I) is 1/DY^2, where DY is the
%           standard deviation associated with Y(I).  If
%           nothing is known about the errors in Y, a con-
%           stant (estimated value) should be used for DY.
%           If PERIOD = TRUE, it is assumed that W(N) =
%           W(1).
%
%       SM = Positive parameter specifying an upper bound on
%            Q2(YS).  H(x) is linear (and Q2 is minimized)
%            if SM is sufficiently large that the constraint
%            is not active.  It is recommended that SM sat-
%            isfy N-SQRT(2N) <= SM <= N+SQRT(2N) and
%            SM = N is reasonable if W(I) = 1/DY^2.
%
%       SMTOL = Parameter in the range (0,1) specifying the
%               relative error allowed in satisfying the
%               constraint:  the constraint is assumed to
%               be satisfied if SM*(1-SMTOL) <= Q2 <=
%               SM*(1+SMTOL).  A reasonable value for SMTOL
%               is SQRT(2/N) for N > 2.
%
% On output:
%
%       YS = Vector of size(X) containing values of H at the
%            abscissae unless IER < 0.  YS(N) = YS(1) if
%            PERIOD = TRUE.
%
%       YP = Vector of size(X) containing first derivative
%            values of H at the abscissae unless IER < 0.
%            YP(N) = YP(1) if PERIOD = TRUE.
%
%       IER = Error indicator:
%             IER = 0 if no errors were encountered and the
%                     constraint is active:  Q2(YS) is ap-
%                     proximately equal to SM.
%             IER = 1 if no errors were encountered but the
%                     constraint is not active:  YS and YP
%                     are the values and derivatives of the
%                     linear function (constant function if
%                     PERIOD = TRUE) that minimizes Q2, and
%                     Q1 = 0.
%             IER = -1 if N, W, SM, or SMTOL is outside its
%                      valid range.  YS and YP are zeros
%                      in this case.
%             IER = -I if X(I) <= X(I-1) for some I in the
%                      range 2 to N.  YS and YP are zeros
%                      in this case.
%
% Modules required by SMCRV:  B2TRI or B2TRIP, SNHCSH,
%                               YPCOEF
%
%***********************************************************

% Set fid = 1 to print diagnostic error messages.

fid = -1;

% Test for errors, and compute the components of the system
%   (normal equations) for the weighted least squares linear
%   fit.

m = size(x);
n = length(x);
if (n < 2  ||  (n < 3  &&  period)  ||  sm <= 0  || ...
    smtol <= 0  ||  smtol >= 1.0  ||  any(w <= 0))
   ys = zeros(m);
   yp = zeros(m);
   ier = -1;
   return;
end

if (any(x(2:n) <= x(1:n-1)))
   ys = zeros(m);
   yp = zeros(m);
   ier = -(find(x(2:n) <= x(1:n-1), 1) + 1);
   return;
end

if (~period)
   t = w.*x;
   c11 = sum(t.*x);   % Sum(w(i)*x(i)^2)
   c12 = sum(t);      % Sum(w(i)*x(i))
   r1 = sum(t.*y);    % Sum(w(i)*x(i)*y(i))
end
c22 = sum(w);         % Sum(w(i))
r2 = sum(w.*y);       % Sum(w(i)*y(i))

% Solve the system for (HP,H0), where HP is the derivative
%   (constant) and H0 = H(0).

if (period)
   h0 = r2/c22;
   hp = 0;
else
   h0 = (c11*r2-c12*r1)/(c11*c22-c12*c12);
   hp = (r1 - c12*h0)/c11;
end

% Store function values and derivatives, and accumulate
%   Q2 = (Y-YS)^T*W*(Y-YS).

ys = hp*x + h0;
yp = hp*ones(m);
t = y-ys;
q2 = dot(t, w.*t);

% Test for the constraint satisfied by the linear fit.

if q2 <= sm*(1.0 + smtol);

%   The constraint is satisfied by a linear function.

   if (fid > 0)
      fprintf(fid,['\n\n\n SMCRV:  The constraint is not ', ...
              'active, and the fit is linear.\n\n']);
   end
   ier = 1;
   return;
end

% Compute the matrix components for the linear system.

ier = 0;
dx = x(2:n)-x(1:n-1);
sig = abs(sigma);
if size(sig) ~= size(dx)
   sig = sig';
end
[d,sd] = ypcoef(sig,dx);

% Compute G0 = G(0), and print a heading.

s = 1.0/sqrt(sm);
g0 = 1.0/sqrt(q2) - s;
if (fid > 0) 
   fprintf(fid,['\n\n\n    SMCRV:  SM = %10.4e, SMTOL = %14.8e, ', ...
           'G(0) = %15.8e\n\n\n'], sm, smtol, g0);
end

% G(P) is strictly increasing and concave, and G(0) < 0.
%
% Initialize parameters for the zero finder.

iter = 0;
tol = min([(1.0-1.0/sqrt(1.0+smtol))*s, (-1.0+1.0/sqrt(1.0-smtol))*s]);
options = optimset('TolX',tol);

% Find a bracketing interval

g = 0;
p = sm;
while g <= 0
   p = 10.0*p;
   g = fsmcrv(p);
end

p = fzero(@fsmcrv,[0,p],options);
if (fid > 0) 
   fprintf(fid,'\n P = %15.8e\n', p);
end
return;

   function g = fsmcrv(p)
% Nested function for evaluation of g.

   if p == 0
      g = g0;
      return;
   end

   if (~period) 
      [ys,yp] = b2tri(x,y,w,p,d,sd);
   else
      [ys,yp] = b2trip(x,y,w,p,d,sd);
   end
   t = y-ys;
   q2 = dot(t, w.*t);
   g = 1.0/sqrt(q2) - s;
   iter = iter + 1;
   if (fid > 0) 
      fprintf(fid,'\n %0.0f:  P = %15.8e,  G = %15.8e\n', ...
              iter, p, g);
   end

   return;
   end
   
end  % smcrv
