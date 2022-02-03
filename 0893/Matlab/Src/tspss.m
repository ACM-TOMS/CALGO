function [sigma,ys,yp,nit,ier,dys,dsmax] = tspss(x,y,per, ...
         sig1,w,sm,smtol)
% tspss:  Parameters defining shape-preserving smoothing curve
%
% USAGE:  [sigma,ys,yp,nit,ier,dys,dsmax] = tspss(x,y,per, ...
%         sig1,w,sm,smtol);
%
%   This function computes a set of parameter values that
% define a smoothing tension spline H(x).  The parameters
% consist of knot values YS and derivatives YP computed
% by Function SMCRV, and tension factors SIGMA computed by
% Function SIGS (unless SIG1 >= 0, indicating uniform 
% tension).  The Hermite interpolatory tension spline H(x) 
% defined by the knot values and derivatives has two contin-
% uous derivatives and satisfies either natural or periodic 
% end conditions.
%
%   The tension spline may be evaluated by Function TSVAL1
% or Functions HVAL (values), HPVAL (first derivatives),
% HPPVAL (second derivatives), HPPPVAL (third derivatives),
% and TSINTL (integrals).
%
% On input:
%
%       X = Vector of length N containing a strictly in-
%           creasing sequence of abscissae:  X(I) < X(I+1)
%           for I = 1 to N-1.  N >= 2 and N >= 3 if PER =
%           TRUE.
%
%       Y = Vector of length N containing data values asso-
%           ciated with the abscissae.  If PER = TRUE, it is
%           assumed that Y(N) = Y(1).
%
%       PER = Logical variable with value TRUE if and only
%             H(x) is to be a periodic function with period
%             X(N)-X(1).  It is assumed without a test that
%             Y(N) = Y(1) in this case.  On output, YP(N) =
%             YP(1) and, more generally, the values and
%             first two derivatives of H at X(1) agree with
%             those at X(N).  If H(x) is one of the compo-
%             nents of a parametric curve, this option may
%             be used to obtained a closed curve.  If PER =
%             FALSE, H satisfies natural end conditions:
%             zero second derivatives at X(1) and X(N).
%
%       SIG1 = Constant (uniform) tension factor in the
%              range 0 to SBIG, or negative value if 
%              variable tension is to be used.  If SIG1 = 0,
%              H(x) is a cubic spline, and as SIG1 
%              increases, H(x) approaches piecewise linear. 
%              If SIG1 < 0, tension factors are chosen (by
%              SIGS) to preserve local monotonicity and
%              convexity of the data.  This may result in a
%              better fit than the case of uniform tension, 
%              but requires an iteration on calls to SMCRV 
%              and SIGS.
%
%       W = Vector of length N containing positive weights
%           associated with the data values.  The recommend-
%           ed value of W(I) is 1/DY^2, where DY is the
%           standard deviation associated with Y(I).  If
%           nothing is known about the errors in Y, a con-
%           stant (estimated value) should be used for DY.
%           If PER = TRUE, it is assumed that W(N) = W(1).
%
%       SM = Positive parameter specifying an upper bound on
%            Q2(YS), where Q2(YS) is the weighted sum of
%            squares of deviations from the data (differ-
%            ences between YS and Y).  H(x) is linear (and
%            Q2 is minimized) if SM is sufficiently large
%            that the constraint is not active.  It is
%            recommended that SM satisfy N-SQRT(2N) <= SM
%            <= N+SQRT(2N) and SM = N is reasonable if
%            W(I) = 1/DY^2.
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
%       SIGMA = Array of size 1 by n-1 containing tension 
%               factors.  SIGMA(I) is associated with inter-
%               val (X(I),X(I+1)) for I = 1 to N-1.  SIGMA 
%               is zeros if N is invalid or -4 < IER < -1,
%               and SIGMA is constant if IER = -1 (and N
%               is valid) or IER = -4.
%
%       YS = Vector of size(X) containing values of H at the
%            abscissae.  YS(N) = YS(1) if PER = TRUE.  YS is
%            zeros if IER < 0.
%
%       YP = Vector of size(X) containing first derivative
%            values of H at the abscissae.  YP(N) = YP(1)
%            if PER = TRUE.  YP is zeros if IER < 0.
%
%       NIT = Number of iterations (calls to SIGS).  NIT = 0
%             if IER < 0 or SIG1 >= 0.  If NIT > 0, NIT+1
%             calls to SMCRV were employed.
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
%                     Q1 = 0 (refer to SMCRV).
%             IER = -1 if N, W, SM, or SMTOL is outside its
%                      valid range.
%             IER = -2 if the number of input arguments is
%                      not valid.
%             IER = -3 if SIG1 > SBIG.
%             IER = -4 if the abscissae X are not strictly
%                      increasing.
%
%       DYS = Maximum relative change in a component of YS
%             on the last iteration if NIT > 0.
%
%       DSMAX = Maximum relative change in a component of
%               SIGMA on the last iteration if NIT > 0.
%
% Modules required by TSPSS:  B2TRI or B2TRIP, SIGS, SMCRV,
%                               SNHCSH, YPCOEF
%
%***********************************************************

global SBIG

stol = 0;
maxit = 99;
dystol = 0.01;

% Convergence parameters:
%
%   STOL = Absolute tolerance for SIGS.
%
%   MAXIT = Maximum number of SMCRV/SIGS iterations.
%
%   DYSTOL = Bound on the maximum relative change in a
%            component of YS defining convergence of
%            the SMCRV/SIGS iteration when SIG1 >= 0.
%

m = size(x);
n = length(x);

sigma = zeros(1,n-1);
nit = 0;
dys = 0;
dsmax = 0;

% Initialize NIT, and test for invalid input parameters N or
%   SIG1.

if (n < 2  ||  (per  &&  n < 3))
   ys = zeros(m);
   yp = zeros(m);
   ier = -1;
   return;
end

% Test for incorrect number of input arguments.

if (nargin ~= 7)
   ys = zeros(m);
   yp = zeros(m);
   ier = -2;
   return;
end

unifrm = (sig1 >= 0);
if (unifrm) 
   sig = sig1;
   if (sig > SBIG)
      ys = zeros(m);
      yp = zeros(m);
      ier = -3;
      return;
   end
else
   sig = 0;
end

% Store uniform tension factors, or initialize SIGMA to
%   zeros.

sigma = sig*ones(1,n-1);

% Compute smoothing curve for uniform tension.

[ys,yp,ier] = smcrv(x,y,sigma,per,w,sm,smtol);
if (ier <= -2), ier = -4; end
if (ier < 0  ||  unifrm), return, end

%   Iterate on calls to SIGS and SMCRV.  The function
%     values YS from the previous iteration are stored 
%     in WK.
%
%   DYS is the maximum relative change in a component of YS.
%   ICNT is the number of tension factors that were
%        increased by SIGS.
%   DSMAX is the maximum relative change in a component of
%         SIGMA.

wk = zeros(m);
e = zeros(m);
i = 2:n-1;
for iter = 1:maxit
   wk(i) = ys(i);
   [sigma,dsmax,icnt] = sigs(x,y,yp,stol,sigma);
   [ys,yp,ierr] = smcrv(x,y,sigma,per,w,sm,smtol);
   e(i) = abs(ys(i)-wk(i));
   k = find(wk);
   e(k) = e(k)/wk(k);
   dys = max(e(i));
   if (icnt == 0  ||  dys <= dystol), break, end
end  % for

% No error encountered.

nit = iter;
ier = ierr;
return;

end  % tspss
