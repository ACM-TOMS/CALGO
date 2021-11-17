function [xp,yp,sigma,ier,dyp,dsmax] = tspbp(t,x,y,ncd, ...
         iendc,per,bl,bu,bmax,xp1,xpn,yp1,ypn)
% tspbp:  Parameters defining constrained parametric planar curve
%
% USAGE:  [xp,yp,sigma,ier,dyp,dsmax] = tspbp(t,x,y,ncd, ...
%         iendc,per,bl,bu,bmax,xp1,xpn,yp1,ypn);
%
%   This function computes a set of values that define a
% parametric planar curve C(t) = (H1(t),H2(t)) whose compo-
% nents are Hermite interpolatory tension splines.  The
% output values consist of knot derivative values XP and YP
% computed by Function YPC1T, YPC2, or YPC2P, and tension
% factors SIGMA chosen (by Function SIGBP) to satisfy user-
% specified bounds on the signed distance between C and the
% polygonal curve associated with the control points (refer
% to BL and BU below).
%
%   Refer to Function TSPSP for an alternative method of
% computing tension factors.
%
%   The tension splines may be evaluated by Function
% TSVAL2 or Functions HVAL (values), HPVAL (first deriva-
% tives), HPPVAL (second derivatives), HPPPVAL (third
% derivatives), and TSINTL (integrals).
%
% On input:
%
%       T = Vector of length N containing a strictly in-
%           creasing sequence of knots (discrete parameter
%           values).  Refer to Function ARCL2D.  N >= 2, and
%           N >= 3 if PER = TRUE.
%
%       X,Y = Vectors of length N containing the Cartesian
%             coordinates of an ordered sequence of data
%             points C(I), I = 1 to N, such that C(I) ~=
%             C(I+1).  C(t) is constrained to pass through
%             these points.  In the case of a closed curve
%             (PER = TRUE), the first and last points should
%             coincide. 
%
%       NCD = Number of continuous derivatives at the knots.
%             NCD = 1 or NCD = 2.  If NCD = 1, XP and YP are
%             the components of unit tangent vectors com-
%             puted as weighted averages of incident chord 
%             directions.  Otherwise, a linear system is 
%             solved for the derivative values required for
%             second derivative continuity.  This requires 
%             iterating on calls to YPC2 or YPC2P and calls 
%             to SIGBP, and generally results in more 
%             nonzero tension factors (hence more expensive 
%             evaluation).
%
%       IENDC = End condition indicator for NCD = 2 and PER
%               = FALSE (or dummy parameter otherwise):
%               IENDC = 0 if XP(1), XP(N), YP(1), and YP(N)
%                         are to be computed by monotonicity-
%                         constrained parabolic fits (YPC1).
%               IENDC = 1 if the first derivatives of H1 at
%                         the left and right endpoints are
%                         user-specified in XP1 and XPN,
%                         respectively, and the first deriv-
%                         atives of H2 at the ends are
%                         specified in YP1 and YPN.
%               IENDC = 2 if the second derivatives of H1
%                         and H2 at the endpoints are user-
%                         specified in XP1, XPN, YP1, and
%                         YPN.
%               IENDC = 3 if the end conditions are to be
%                         computed by Function ENDSLP and
%                         vary with SIGMA(1) and SIGMA(N-1).
%
%       PER = Logical variable with value TRUE if and only
%             a closed curve is to be constructed:  H1(t)
%             and H2(t) are to be periodic functions with
%             period T(N)-T(1), where T(1) and T(N) are the
%             parameter values associated with the first and
%             last data points.  It is assumed that X(N) =
%             X(1) and Y(N) = Y(1) in this case, and, on
%             output, XP(N) = XP(1) and YP(N) = YP(1).
%
%       BL,BU = Vectors of length N-1 containing (for each
%               knot subinterval [t1,t2]) lower and upper 
%               bounds, respectively, on the signed perpen-
%               dicular distance d(t) = (C2-C1)/DC X (C(t)-
%               C1), where C1 and C2 are the ordered data 
%               points associated with the interval, and DC
%               is the length of the line segment C1-C2,
%               assumed equal to t2-t1.  If the curve is not
%               parameterized by cumulative arc length
%               (Function ARCL2D), then BL and BU should be
%               input as the required bounds scaled by DC/
%               (t2-t1).  Note that d(t) > 0 iff C(t) lies
%               strictly to the left of the line segment as
%               viewed from C1 toward C2.  For I = 1 to N-1,
%               SIGMA(I) is chosen to be as small as possible
%               within the constraint that BL(I) <= d(t) <= 
%               BU(I) for all t in the interval.  BL(I) < 0 
%               and BU(I) > 0 for I = 1 to N-1.  A null 
%               constraint is specified by BL(I) <= -BMAX or 
%               BU(I) >= BMAX.
%
%       BMAX = User-defined value of infinity which, when
%              used as an upper bound in BU (or when its
%              negative is used as a lower bound in BL),
%              specifies that no constraint is to be en-
%              forced.
%
%       XP1,XPN = End condition values if NCD = 2 and
%                 IENDC = 1 or IENDC = 2.
%
%       YP1,YPN = End condition values if NCD = 2 and
%                 IENDC = 1 or IENDC = 2.
%
% On output:
%
%       XP = Array of size(T) containing derivatives of H1 
%            at the knots.  XP is zeros if -4 < IER < 0,
%            and XP is only partially computed if IER = -5.
%
%       YP = Array of size(T) containing derivatives of H2 
%            at the knots.  YP is zeros if -4 < IER < 0,
%            and YP is only partially computed if IER = -5.
%
%       SIGMA = Array of size 1 by N-1 containing tension 
%               factors for which C(t) satisfies the con-
%               straints defined by BL and BU.  SIGMA(I) is
%               associated with interval (T(I),T(I+1)) for 
%               I = 1 to N-1.  SIGMA(I) is limited to SBIG 
%               (in which case C(t) is close to the line 
%               segment associated with the interval), and 
%               if no constraint is specified in the 
%               interval, then SIGMA(I) = 0, and thus
%               H1 and H2 are cubic functions of t.  SIGMA 
%               is zeros if IER < 0 and IER ~= -3.
%
%       IER = Error indicator or iteration count:
%             IER = IC >= 0 if no errors were encountered
%                      and IC calls to SIGBP and IC+1 calls
%                      to YPC1, YPC1P, YPC2 or YPC2P were
%                      employed.  (IC = 0 if NCD = 1).
%             IER = -1 if N, NCD, or IENDC is outside its
%                      valid range.
%             IER = -2 if the number of input arguments is
%                      not consistent with the values.
%             IER = -3 if BL(I) >= 0 or BU(I) <= 0 for
%                      some I in the range 1 to N-1.
%                      SIGMA(J) = 0 for J >= I in this
%                      case.
%             IER = -4 if a pair of adjacent data points
%                      coincide:  X(I) = X(I+1) and Y(I) =
%                      Y(I+1) for some I in the range 1 to
%                      N-1.
%             IER = -5 if the knots are not strictly 
%                      increasing.
%
%       DYP = Maximum relative change in a component of XP 
%             or YP on the last iteration if IER > 0.
%
%       DSMAX = Maximum relative change in a component of
%               SIGMA on the last iteration if IER > 0.
%
% Modules required by TSPBP:  ENDSLP, SIGBP, SNHCSH, YPCOEF,
%                               YPC1T, YPC2, YPC2P
%
%***********************************************************

stol = 0;
maxit = 49;
dyptol = 0.01;

% Convergence parameters:
%
%   STOL = Absolute tolerance for SIGBP.
%
%   MAXIT = Maximum number of YPC2/SIGBP iterations for each
%           loop if NCD = 2.
%
%   DYPTOL = Bound on the maximum relative change in a
%            component of XP or YP defining convergence
%            of the YPC2/SIGBP iteration when NCD = 2.

m = size(t);
n = length(t);

sigma = zeros(1,n-1);
dyp = 0;
dsmax = 0;

% Test for invalid input parameters N or NCD.

if (n < 2  ||  (per  &&  n < 3)  ||  ncd < 1  ||  ncd > 2)
   xp = zeros(m);
   yp = zeros(m);
   ier = -1;
   return;
end

% Test for incorrect number of input arguments.

userc = ~per  &&  ncd == 2  &&  (iendc == 1  ||  iendc == 2);
if (~userc  &&  nargin ~= 9)  ||  (userc  &&  nargin ~= 13)
   xp = zeros(m);
   yp = zeros(m);
   ier = -2;
   return;
end

% Initialize iteration count ITER.

iter = 0;
if (ncd == 1)

% NCD = 1.

   [xp,yp,ierr] = ypc1t(x,y);
   if (ierr ~= 0) 
      ier = -4;
      return;
   end

%   Compute tension factors.

   [sigma,dsmax,ierr] = sigbp(x,y,xp,yp,stol,bl,bu,bmax,sigma);
   if (ierr < 0) 
      ier = -3;
   else
      ier = 0;
   end
   return;
end

% NCD = 2.

if (~per) 

%   Nonperiodic case:  call YPC2 and test for IENDC invalid.

   if (iendc == 1  ||  iendc == 2)
      [xp,ierx] = ypc2(t,x,sigma,iendc,iendc,xp1,xpn);
      [yp,iery] = ypc2(t,y,sigma,iendc,iendc,yp1,ypn);
   else
      [xp,ierx] = ypc2(t,x,sigma,iendc,iendc);
      [yp,iery] = ypc2(t,y,sigma,iendc,iendc);
   end
   if (ierx == 1  ||  iery == 1) 
      ier = -1;
      return;
   end
   if (ierx > 1  ||  iery > 1)
      ier = -5;
      return;
   end
else

%   Periodic fit:  call YPC2P.

   [xp,ierx] = ypc2p(t,x,sigma);
   [yp,iery] = ypc2p(t,y,sigma);
   if (ierx ~= 0  ||  iery ~= 0) 
      ier = -5;
      return;
   end
end

%   Iterate on calls to SIGBP and YPC2 (or YPC2P).  The
%     derivative estimates XP and YP from the previous 
%     iteration are stored in WX and WY.
%
%   LOOP2 is TRUE iff tension factors are not allowed to
%         decrease between iterations (loop 1 failed to
%         converge with MAXIT iterations).
%   DYP is the maximum relative change in a component of XP
%       or YP.
%   ICNT is the number of tension factors that were altered
%        by SIGBP.
%   DSMAX is the maximum relative change in a component of
%         SIGMA.

loop2 = false;
while (true)
   for iter = 1:maxit
      i = 2:n-1;
      wx(i) = xp(i);
      wy(i) = yp(i);
      [sigma,dsmax,icnt] = sigbp(x,y,xp,yp,stol,bl,bu,bmax,sigma);
      if (icnt < 0) 
         ier = -3;
         return;
      end
      if (~per) 
         if (iendc == 1  ||  iendc == 2)
            xp = ypc2(t,x,sigma,iendc,iendc,xp1,xpn);
            yp = ypc2(t,y,sigma,iendc,iendc,yp1,ypn);
         else
            xp = ypc2(t,x,sigma,iendc,iendc);
            yp = ypc2(t,y,sigma,iendc,iendc);
         end
      else
         xp = ypc2p(t,x,sigma);
         yp = ypc2p(t,y,sigma);
      end
      ex(i) = abs(xp(i)-wx(i));
      k = find(wx);
      ex(k) = ex(k)/abs(wx(k));
      ey(i) = abs(yp(i)-wy(i));
      k = find(wy);
      ey(k) = ey(k)/abs(wy(k));
      dyp = max([ex(i), ey(i)]);
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
    
end  % tspbp
