function [yp,sigma,icflg,ier,dyp,dsmax] = tspbi(x,y,ncd, ...
         iendc,per,b,bmax,yp1,ypn)
% tspbi:  Parameters defining constrained interpolatory tension spline
%
% USAGE:  [yp,sigma,icflg,ier,dyp,dsmax] = tspbi(x,y,ncd, ...
%         iendc,per,b,bmax,yp1,ypn);
%
%   This function computes a set of parameter values that
% define a Hermite interpolatory tension spline H(x).  The
% parameters consist of knot derivative values YP computed
% by Function YPC1, YPC1P, YPC2, or YPC2P, and tension
% factors SIGMA chosen to satisfy user-specified constraints
% (by Function SIGBI).  Refer to Function TSPSI for an
% alternative method of computing tension factors.
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
%           I = 1 to N.  If NCD = 1 and PER = TRUE, Y(1)
%           and Y(N) should be identical.
%
%       NCD = Number of continuous derivatives at the knots.
%             NCD = 1 or NCD = 2.  If NCD = 1, the YP values
%             are computed by local monotonicity-constrained
%             quadratic fits.  Otherwise, a linear system is
%             solved for the derivative values that result
%             in second derivative continuity.  This re-
%             quires iterating on calls to YPC2 or YPC2P and
%             calls to SIGBI, and generally results in more
%             nonzero tension factors (hence more expensive
%             evaluation).
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
%       B = Array dimensioned 5 by N-1 containing bounds or
%           flags which define the constraints.  For I = 1
%           to N-1, column I defines the constraints associ-
%           ated with interval (X(I),X(I+1)) as follows:
%
%             B(1,I) is an upper bound on H
%             B(2,I) is a lower bound on H
%             B(3,I) is an upper bound on HP
%             B(4,I) is a lower bound on HP
%             B(5,I) specifies the required sign of HPP
%
%           where HP and HPP denote the first and second
%           derivatives of H, respectively.  A null con-
%           straint is specified by abs(B(K,I)) >= BMAX
%           for K < 5, or B(5,I) = 0:   B(1,I) >= BMAX,
%           B(2,I) <= -BMAX, B(3,I) >= BMAX, B(4,I) <=
%           -BMAX, or B(5,I) = 0.  Any positive value of
%           B(5,I) specifies that H should be convex, a
%           negative values specifies that H should be con-
%           cave, and 0 specifies that no restriction be
%           placed on HPP.  Refer to Functions SIG0, SIG1,
%           and SIG2 for definitions of valid constraints.
%
%       BMAX = User-defined value of infinity which, when
%              used as an upper bound in B (or when 
%              its negative is used as a lower bound),
%              specifies that no constraint is to be en-
%              forced.
%
%       YP1,YPN = End condition values if NCD = 2 and
%                 IENDC = 1 or IENDC = 2.
%
% On output:
%
%       YP = Array of size(X) containing derivatives of H at 
%            the abscissae.  YP is zeros if -3 < IER < 0,
%            and YP is only partially computed if IER = -4.
%
%       SIGMA = Array of size 1 by N-1 containing tension 
%               factors for which H(x) satisfies the 
%               constraints defined by B.  SIGMA(I) is 
%               associated with interval (X(I),X(I+1)) for 
%               I = 1 to N-1.  If infinite tension is 
%               required in interval I, then SIGMA(I) = SBIG
%               (and H is an approximation to the linear 
%               interpolant on the interval), and if no 
%               constraint is specified in the interval, 
%               then SIGMA(I) = 0, and thus H is cubic.  
%               Invalid constraints are treated as null 
%               constraints.  SIGMA is zeros if IER < 0.
%
%       ICFLG = Array of size 1 by N-1 containing invalid 
%               constraint flags associated with intervals.
%               For I = 1 to N-1, ICFLG(I) is a 5-bit value
%               b5b4b3b2b1, where bK = 1 if and only if 
%               constraint K cannot be satisfied.  Thus, all
%               constraints in interval I are satisfied if 
%               and only if ICFLG(I) = 0 (and IER >= 0).
%               ICFLG is zeros if IER < 0.
%
%       IER = Error indicator or iteration count:
%             IER = IC >= 0 if no errors were encountered
%                      (other than invalid constraints) and
%                      IC calls to SIGBI and IC+1 calls to
%                      YPC1, YPC1P, YPC2 or YPC2P were
%                      employed.  (IC = 0 if NCD = 1).
%             IER = -1 if N, NCD, or IENDC is outside its
%                      valid range.
%             IER = -2 if the number of input arguments is
%                      not consistent with the values.
%             IER = -4 if the abscissae X are not strictly
%                      increasing.
%
%       DYP = Maximum relative change in a component of YP 
%             on the last iteration if IER > 0.
%
%       DSMAX = Maximum relative change in a component of
%               SIGMA on the last iteration if IER > 0.
%
% Modules required by TSPBI:  ENDSLP, SIG0, SIG1, SIG2,
%                               SIGBI, SNHCSH, YPCOEF, YPC1,
%                               YPC1P, YPC2, YPC2P
%
%***********************************************************

stol = 0;  
maxit = 49;  
dyptol = 0.01;

% Convergence parameters:
%
%   STOL = Absolute tolerance for SIGBI.
%
%   MAXIT = Maximum number of YPC2/SIGBI iterations for
%           each loop if NCD = 2.
%
%   DYPTOL = Bound on the maximum relative change in a
%            component of YP defining convergence of
%            the YPC2/SIGBI iteration when NCD = 2.

m = size(x);
n = length(x);

sigma = zeros(1,n-1);
icflg = zeros(1,n-1);
dyp = 0;
dsmax = 0;

% Test for invalid input parameters N or NCD.

if (n < 2  ||  (per  &&  n < 3)  ||  ncd < 1  ||  ncd > 2) 
   yp = zeros(m);
   ier = -1;
   return;
end

% Test for incorrect number of input arguments.

userc = ~per  &&  ncd == 2  &&  (iendc == 1  ||  iendc == 2);
if (~userc  &&  nargin ~= 7)  ||  (userc  &&  nargin ~= 9)
   yp = zeros(m);
   ier = -2;
   return;
end

% Initialize iteration count ITER.

iter = 0;
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

%   Compute tension factors.

   [sigma,icflg,dsmax] = sigbi(x,y,yp,stol,b,bmax,sigma);
   ier = iter;
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
      return
   end
end

%   Iterate on calls to SIGBI and YPC2 (or YPC2P).  The
%     derivative estimates YP from the previous iteration
%     are stored in WK.
%
%   LOOP2 is TRUE iff tension factors are not allowed to
%         decrease between iterations (loop 1 failed to
%         converge with MAXIT iterations).
%   DYP is the maximum relative change in a component of YP.
%   ICNT is the number of tension factors that were altered
%        by SIGBI.
%   DSMAX is the maximum relative change in a component of
%         SIGMA.

wk = zeros(m);
e = zeros(m);
i = 2:n-1;
loop2 = false;
while (true)
   for iter = 1:maxit
      wk(i) = yp(i);
      [sigma,icflg,dsmax,icnt] = sigbi(x,y,yp,stol,b,bmax,sigma);
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
      if (~loop2) 

%   Loop 1:  reinitialize SIGMA to zeros.

         sigma = zeros(1,n-1);
      end
   end  % for

%   The loop failed to converge within MAXIT iterations.

   if (loop2), break, end
   loop2 = true;
end  % while

% Update iter.

if (loop2), iter = iter + maxit; end

% No error encountered.

ier = iter;
return;

end  % tspbi
