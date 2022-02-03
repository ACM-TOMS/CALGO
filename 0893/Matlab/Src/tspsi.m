function [yp,sigma,ier,dyp,dsmax] = tspsi(x,y,ncd,iendc, ...
         per,sig1,yp1,ypn)
% tspsi:  Parameters defining shape-preserving interpolatory tension spline
%
% USAGE:  [yp,sigma,ier,dyp,dsmax] = tspsi(x,y,ncd,iendc, ...
%         per,sig1,yp1,ypn);
%
%
%   This function computes a set of parameter values that
% define a Hermite interpolatory tension spline H(x).  The
% parameters consist of knot derivative values YP computed
% by Function YPC1, YPC1P, YPC2, or YPC2P, and tension
% factors SIGMA computed by Function SIGS (unless SIG1 >= 
% 0, indicating uniform tension).  Alternative methods for 
% computing SIGMA are provided by Function TSPBI and Func-
% tions SIG0, SIG1, and SIG2.
%
%   Refer to Function TSPSS for a means of computing
% parameters that define a smoothing curve rather than an
% interpolatory curve.
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
%           ciated with the abscissae.  H(X(I)) = Y(I) for
%           I = 1 to N.  If NCD = 1 and PER = TRUE, Y(1) and
%           Y(N) should be identical.
%
%       NCD = Number of continuous derivatives at the knots.
%             NCD = 1 or NCD = 2.  If NCD = 1, the YP values
%             are computed by local monotonicity-constrained
%             quadratic fits.  Otherwise, a linear system is
%             solved for the derivative values that result
%             in second derivative continuity.  Unless
%             SIG1 >= 0, this requires iterating on calls
%             to YPC2 or YPC2P and calls to SIGS, and
%             generally results in more nonzero tension
%             factors (hence more expensive evaluation).
%
%       IENDC = End condition indicator for NCD = 2 and PER
%               = FALSE (or dummy parameter otherwise):
%               IENDC = 0 if YP(1) and YP(N) are to be com-
%                         puted by monotonicity-constrained
%                         parabolic fits to the first three
%                         and last three points, respective-
%                         ly.  This is identical to the
%                         values computed by YPC1.
%               IENDC = 1 if the first derivatives of H at
%                         X(1) and X(N) are user-specified
%                         in YP1 and YPN, respectively.
%               IENDC = 2 if the second derivatives of H at
%                         X(1) and X(N) are user-specified
%                         in YP1 and YPN, respectively.
%               IENDC = 3 if the end conditions are to be
%                         computed by Function ENDSLP and
%                         vary with SIGMA(1) and SIGMA(N-1).
%
%       PER = Logical variable with value TRUE if and only
%             H(x) is to be a periodic function with period
%             X(N)-X(1).  It is assumed without a test that
%             Y(N) = Y(1) in this case.  On output, YP(N) =
%             YP(1).  If H(x) is one of the components of a
%             parametric curve, this option may be used to
%             obtained a closed curve.
%
%       SIG1 = Constant (uniform) tension factor in the
%              range 0 to SBIG, or negative value if 
%              variable tension is to be used.  If SIG1 = 0,
%              H(x) is piecewise cubic (a cubic spline if
%              NCD = 2), and as SIG1 increases, H(x) 
%              approaches the piecewise linear interpolant.
%              If SIG1 < 0, tension factors are chosen (by
%              SIGS) to preserve local monotonicity and
%              convexity of the data.  This often improves
%              the appearance of the curve over the piece-
%              wise cubic fit.
%
%       YP1,YPN = End condition values if NCD = 2 and
%                 IENDC = 1 or IENDC = 2.
%
% On output:
%
%       YP = Array of size(X) containing derivatives of H at
%            the abscissae.  YP is zeros if -4 < IER < 0,
%            and YP is only partially computed if IER = -4.
%
%       SIGMA = Array of size 1 by N-1 containing tension 
%               factors.  SIGMA(I) is associated with inter-
%               val (X(I),X(I+1)) for I = 1 to N-1.  SIGMA 
%               is zeros if -4 < IER < 0 (unless IENDC
%               is invalid), and SIGMA is constant (not 
%               optimal) if IER = -4 or IENDC (if used) is 
%               invalid.
%
%       IER = Error indicator or iteration count:
%             IER = IC >= 0 if no errors were encountered
%                      and IC calls to SIGS and IC+1 calls
%                      to YPC1, YPC1P, YPC2 or YPC2P were
%                      employed.  (IC = 0 if NCD = 1).
%             IER = -1 if N, NCD, or IENDC is outside its
%                      valid range.
%             IER = -2 if the number of input arguments is
%                      not consistent with the values.
%             IER = -3 if SIG1 > SBIG.
%             IER = -4 if the abscissae X are not strictly
%                      increasing.
%
%       DYP = Maximum relative change in a component of YP 
%             on the last iteration if IER > 0.
%
%       DSMAX = Maximum relative change in a component of
%               SIGMA on the last iteration if IER > 0.
%
% Modules required by TSPSI:  ENDSLP, SIGS, SNHCSH, YPCOEF,
%                               YPC1, YPC1P, YPC2, YPC2P
%
%***********************************************************

global SBIG

stol = 0;
maxit = 99;
dyptol = 0.01;

% Convergence parameters:
%
%   STOL = Absolute tolerance for SIGS.
%
%   MAXIT = Maximum number of YPC2/SIGS iterations.
%
%   DYPTOL = Bound on the maximum relative change in a
%            component of YP defining convergence of
%            the YPC2/SIGS iteration when NCD = 2 and
%            SIG1 < 0.

m = size(x);
n = length(x);

dyp = 0;
dsmax = 0;

% Test for invalid input parameters N and NCD.

if (n < 2  ||  (per  &&  n < 3)  ||  ncd < 1  ||  ncd > 2)
   yp = zeros(m);
   sigma = zeros(1,n-1);
   ier = -1;
   return;
end

% Test for incorrect number of input arguments.

userc = ~per  &&  ncd == 2  &&  (iendc == 1  ||  iendc == 2);
if (~userc  &&  nargin ~= 6)  ||  (userc  &&  nargin ~= 8)
   yp = zeros(m);
   sigma = zeros(1,n-1);
   ier = -2;
   return;
end

if (sig1 >= 0) 
   sig = sig1;
   if (sig > SBIG)
      yp = zeros(m);
      sigma = zeros(1,n-1);
      ier = -3;
      return;
   end
else
   sig = 0;
end

% Initialize iteration count ITER, and store uniform
%   tension factors, or initialize SIGMA to zeros.

iter = 0;
unifrm = sig1 >= 0;
sigma = sig*ones(1,n-1);
if (ncd == 1) 

% NCD = 1.

   if (~per) 
      [yp,ierr] = ypc1(x,y);
   else
      [yp,ierr] = ypc1p(x,y);
   end
   if (ierr ~= 0) 
      ier = -4;
      return;
   end
   if (~unifrm)

%   Call SIGS for SIG1 < 0.

      [sigma,dsmax] = sigs(x,y,yp,stol,sigma);
   end
   ier = 0;
   return;
end

% NCD = 2.

if (~per) 

%   Nonperiodic case:  call YPC2 and test for IENDC or X
%     invalid.

   if (iendc == 1  ||  iendc == 2)
      [yp,ierr] = ypc2(x,y,sigma,iendc,iendc,yp1,ypn);
   else
      [yp,ierr] = ypc2(x,y,sigma,iendc,iendc);
   end
   if (ierr == 1) 
      ier = -1;
      return;
   end
   if (ierr > 1)
      ier = -4;
      return;
   end
else

%   Periodic fit:  call YPC2P.

   [yp,ierr] = ypc2p(x,y,sigma);
   if (ierr > 1) 
      ier = -4;
      return;
   end
end
if (unifrm)
   ier = 0;
   return;
end

%   Iterate on calls to SIGS and YPC2 (or YPC2P).  The 
%     derivative estimates YP from the previous iteration
%     are stored in WK.
%
%   DYP is the maximum relative change in a component of YP.
%   ICNT is the number of tension factors that were
%        increased by SIGS.
%   DSMAX is the maximum relative change in a component of
%         SIGMA.

wk = zeros(m);
e = zeros(m);
i = 2:n-1;
for iter = 1:maxit
   wk(i) = yp(i);
   [sigma,dsmax,icnt] = sigs(x,y,yp,stol,sigma);
   if (~per) 
      if (iendc == 1  ||  iendc == 2)
         yp = ypc2(x,y,sigma,iendc,iendc,yp1,ypn);
      else
         yp = ypc2(x,y,sigma,iendc,iendc);
      end
   else
      yp = ypc2p(x,y,sigma);
   end
   e(i) = abs(yp(i)-wk(i));
   k = find(wk);
   e(k) = e(k)/abs(wk(k));
   dyp = max(e(i));
   if (icnt == 0  ||  dyp <= dyptol), break, end
end  % for

% No error encountered.

ier = iter;
return;

end  % tspsi
