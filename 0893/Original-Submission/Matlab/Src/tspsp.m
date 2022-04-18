function [xp,yp,zp,sigma,ier,dyp,dsmax] = tspsp(nd,t,x,y,z, ...
         ncd,iendc,per,sig1,xp1,xpn,yp1,ypn,zp1,zpn)
% tspsp:  Parameters defining shape-preserving parametric space curve
%
% USAGE:  [xp,yp,zp,sigma,ier,dyp,dsmax] = tspsp(nd,t,x,y,z, ...
%         ncd,iendc,per,sig1,xp1,xpn,yp1,ypn,zp1,zpn);
%
%
%   This function computes a set of values that define a
% parametric planar curve C(t) = (H1(t),H2(t)) or space
% curve C(t) = (H1(t),H2(t),H3(t)) whose components are 
% Hermite interpolatory tension splines.  The output values
% consist of knot derivative values XP, YP, (and ZP) 
% computed by Function YPC1T, YPC2, or YPC2P, and tension 
% factors SIGMA computed by Function SIGSP (unless SIG1 >= 
% 0, indicating uniform tension).
%
%   Refer to Function TSPBP for an alternative method of
% computing tension factors in the case of a planar curve.
%
%   The tension splines may be evaluated by Function
% TSVAL2 (or TSVAL3) or Functions HVAL (values), HPVAL
% (first derivatives), HPPVAL (second derivatives), 
% HPPPVAL (third derivatives), and TSINTL (integrals).
%
% On input:
%
%       ND = Number of dimensions:
%            ND = 2 if a planar curve is to be constructed.
%            ND = 3 if a space curve is to be constructed.
%
%       T = Vector of length N containing a strictly in-
%           creasing sequence of knots (discrete parameter
%           values).  Refer to Function ARCL2D and ARCL3D.
%           N >= 2, and N >= 3 if PER = TRUE.
%
%       X,Y,Z = Vectors of length N containing the Cartesian
%               coordinates of an ordered sequence of data
%               points C(I), I = 1 to N, such that C(I) ~=
%               C(I+1).  The curve is constrained to pass
%               through these points.  Z is an unused dummy 
%               parameter if ND = 2.  If PER = TRUE (closed
%               curve), the first and last points should 
%               coincide.
%
%       NCD = Number of continuous derivatives at the knots.
%             NCD = 1 or NCD = 2.  If NCD = 1, XP, YP, (and 
%             ZP) are the components of unit tangent vectors 
%             computed as weighted averages of incident
%             chord directions.  Otherwise, a linear system
%             is solved for the derivative values required
%             for second derivative continuity.  Unless 
%             SIG1 >= 0, this requires iterating on calls to
%             YPC2 or YPC2P and calls to SIGSP, and generally
%             results in more nonzero tension factors (hence
%             more expensive evaluation).
%
%       IENDC = End condition indicator for NCD = 2 and PER
%               = FALSE (or dummy parameter otherwise):
%               IENDC = 0 if XP(1), XP(N), YP(1), YP(N) (and
%                         ZP(1) and ZP(N)) are to be com-
%                         puted by monotonicity-constrained
%                         parabolic fits (YPC1).
%               IENDC = 1 if the first derivatives of H1 at
%                         the left and right endpoints are
%                         user-specified in XP1 and XPN,
%                         respectively, the first deriva-
%                         tives of H2 at the ends are
%                         specified in YP1 and YPN, and,
%                         if ND = 3, the first derivatives
%                         of H3 are specified in ZP1 and
%                         ZPN.
%               IENDC = 2 if the second derivatives of H1,
%                         H2, (and H3) at the endpoints are
%                         user-specified in XP1, XPN, YP1,
%                         YPN, (ZP1, and ZPN).
%               IENDC = 3 if the end conditions are to be
%                         computed by Function ENDSLP and
%                         vary with SIGMA(1) and SIGMA(N-1).
%
%       PER = Logical variable with value TRUE if and only
%             a closed curve is to be constructed:  H1(t),
%             H2(t), (and H3(t)) are to be periodic func-
%             tions with period T(N)-T(1), where T(1) and
%             T(N) are the parameter values associated with
%             the first and last data points.  It is assumed
%             in this case that X(N) = X(1), Y(N) = Y(1)
%             and, if ND = 3, Z(N) = Z(1), and, on output,
%             XP(N) = XP(1), YP(N) = YP(1), (and ZP(N) =
%             ZP(1) if ND = 3).
%
%       SIG1 = Constant (uniform) tension factor in the
%              range 0 to SBIG, or negative value if 
%              variable tension is to be used.  If SIG1 = 0,
%              H(t) is piecewise cubic (a cubic spline if
%              if NCD = 2), and as SIG1 increases, H(t) 
%              approaches the piecewise linear interpolant, 
%              where H is H1, H2, or H3.  If SIG1 < 0,
%              tension factors are chosen (by SIGSP) to
%              preserve local convexity of the data.
%
%       XP1,XPN = End condition values if NCD = 2 and
%                 IENDC = 1 or IENDC = 2.
%
%       YP1,YPN = End condition values if NCD = 2 and
%                 IENDC = 1 or IENDC = 2.
%
%       ZP1,ZPN = End condition values if ND = 3 and NCD = 2
%                 and IENDC = 1 or IENDC = 2.
%
% On output:
%
%       XP = Array of size(T) containing derivatives of H1 
%            at the knots.  XP is zeros if -4 < IER < 0,
%            and XP is only partially computed if IER = -5.
%
%       YP = Array of size(T) containing derivatives of H2 
%            at the knots.  YP is zeros if -4 < IER  < 0,
%            and YP is only partially computed if IER = -5.
%
%       ZP = Array of size(T) containing derivatives of H3 
%            at the knots if ND = 3.  ZP is zeros if -4 < 
%            IER < 0, and ZP is only partially computed if 
%            IER = -5.
%
%       SIGMA = Array of size 1 by N-1 containing tension 
%               factors.  SIGMA(I) is associated with inter-
%               val (T(I),T(I+1)) for I = 1 to N-1.  SIGMA 
%               is zeros if -3 <= IER < 0 (unless IENDC is
%               invalid), and SIGMA is constant (not
%               optimal) if IENDC (if used) is invalid or
%               IER = -4 or IER = -5.
%
%       IER = Error indicator or iteration count:
%             IER = IC >= 0 if no errors were encountered
%                      and IC calls to SIGSP and IC+1 calls
%                      to YPC2 or YPC2P were employed.  
%                      (IC = 0 if NCD = 1).
%             IER = -1 if ND, N, NCD, or IENDC is outside 
%                      its valid range.
%             IER = -2 if the number of input arguments is
%                      not consistent with the values.
%             IER = -3 if SIG1 > SBIG.
%             IER = -4 if a pair of adjacent control points
%                      coincide.
%             IER = -5 if the knots are not strictly
%                      increasing.
%
%       DYP = Maximum relative change in a component of XP,
%             YP, or ZP on the last iteration if IER > 0.
%
%       DSMAX = Maximum relative change in a component of
%               SIGMA on the last iteration if IER > 0.
%
% Modules required by TSPSP:  ENDSLP, SIGSP, SNHCSH, YPCOEF,
%                               YPC1T, YPC2, YPC2P
%
%***********************************************************

global SBIG

stol = 0;
maxit = 99;
dyptol = 0.01;

% Convergence parameters:
%
%   STOL = Absolute tolerance for SIGSP.
%  
%   MAXIT = Maximum number of YPC2/SIGSP iterations.
%
%   DYPTOL = Bound on the maximum relative change in a com-
%            ponent of XP, YP, or ZP defining convergence
%            of the YPC2/SIGSP iteration when NCD = 2 and
%            SIG1 < 0.

m = size(t);
n = length(t);

sigma = zeros(1,n-1);
zp = zeros(m);
dyp = 0;
dsmax = 0;

% Test for invalid input parameters N, ND, or NCD.

if (n < 2  ||  (per  &&  n < 3)  ||  nd < 2  ||  nd > 3  || ...
    ncd < 1  ||  ncd > 2)
   xp = zeros(m);
   yp = zeros(m);
   ier = -1;
   return;
end

% Test for incorrect number of input arguments.

userc = ~per  &&  ncd == 2  &&  (iendc == 1  ||  iendc == 2);
if (~userc  &&  nargin ~= 9)  ||  (userc  &&  nargin ~= 9+2*nd)
   xp = zeros(m);
   yp = zeros(m);
   ier = -2;
   return;
end

unifrm = (sig1 >= 0);
if (unifrm) 
   sig = sig1;
   if (sig > SBIG) 
      xp = zeros(m);
      yp = zeros(m);
      ier = -3;
      return;
   end
else
   sig = 0;
end

% Set SCURV (TRUE iff space curve), initialize iteration 
% count ITER, and store uniform tension factors SIGMA.

scurv = (nd == 3);
iter = 0;
sigma = sig*ones(1,n-1);
ierz = 0;
if (ncd == 1)

% NCD = 1.

   if scurv
      [xp,yp,zp,ierr] = ypc1t(x,y,z);
   else
      [xp,yp,ierr] = ypc1t(x,y);
   end
   if (ierr ~= 0)
      ier = -4;
      return;
   end
   if (~unifrm) 

%   Call SIGSP for UNIFRM = FALSE (SIG1 < 0).

      [sigma,dsmax,ierr] = sigsp(nd,t,x,y,z,xp,yp,zp,stol,sigma);
      if (ierr < -1) 
         ier = -5;
         return;
      end
   end
   ier = 0;
   return;
end

% NCD = 2.

if (~per) 

%   Nonperiodic case:  call YPC2 and test for IENDC invalid.

   if (iendc == 1  ||  iendc == 2)
      [xp,ierx] = ypc2(t,x,sigma,iendc,iendc,xp1,xpn);
      [yp,iery] = ypc2(t,y,sigma,iendc,iendc,yp1,ypn);
      if (scurv)
         [zp,ierz] = ypc2(t,z,sigma,iendc,iendc,zp1,zpn);
      end
   else
      [xp,ierx] = ypc2(t,x,sigma,iendc,iendc);
      [yp,iery] = ypc2(t,y,sigma,iendc,iendc);
      if (scurv)
         [zp,ierz] = ypc2(t,z,sigma,iendc,iendc);
      end
   end
   if (ierx == 1  ||  iery == 1  ||  ierz == 1)
      ier = -1;
      return;
   end
   if (ierx > 1  ||  iery > 1  ||  ierz > 1)
      ier = -5;
      return;
   end
else

%   Periodic fit:  call YPC2P.

   [xp,ierx] = ypc2p(t,x,sigma);
   [yp,iery] = ypc2p(t,y,sigma);
   if (scurv), [zp,ierz] = ypc2p(t,z,sigma); end
   if (ierx ~= 0  ||  iery ~= 0  ||  ierz ~= 0)
      ier = -5;
      return;
   end
end
if (unifrm)
   ier = 0;
   return;
end

%   Iterate on calls to SIGSP and YPC2 (or YPC2P).  The
%     derivative estimates XP, YP, (and ZP) from the
%     previous iteration are stored in WX, WY, (and WZ).
%
%   DYP is the maximum relative change in a component of XP,
%       YP, or ZP.
%   ICNT is the number of tension factors that were
%        increased by SIGSP.
%   DSMAX is the maximum relative change in a component of
%         SIGMA.

for iter = 1:maxit
   wx = xp;
   wy = yp;
   if (scurv), wz = zp; end
   [sigma,dsmax,icnt] = sigsp(nd,t,x,y,z,xp,yp,zp,stol,sigma);
   if (~per)
      if (iendc == 1  ||  iendc == 2)
         xp = ypc2(t,x,sigma,iendc,iendc,xp1,xpn);
         yp = ypc2(t,y,sigma,iendc,iendc,yp1,ypn);
         if (scurv)
            zp = ypc2(t,z,sigma,iendc,iendc,zp1,zpn);
         end
      else
         xp = ypc2(t,x,sigma,iendc,iendc);
         yp = ypc2(t,y,sigma,iendc,iendc);
         if (scurv)
            zp = ypc2(t,z,sigma,iendc,iendc);
         end
      end
   else
      xp = ypc2p(t,x,sigma);
      yp = ypc2p(t,y,sigma);
      if (scurv), zp = ypc2p(t,z,sigma); end
   end
   ex = abs(xp-wx);
   k = find(wx);
   ex(k) = ex(k)./abs(wx(k));
   ey = abs(yp-wy);
   k = find(wy);
   ey(k) = ey(k)./abs(wy(k));
   if (scurv)
      ez = abs(zp-wz);
      k = find(wz);
      ez(k) = ez(k)./abs(wz(k));
   else
      ez = 0;
   end
   dyp = max([ex(:); ey(:); ez(:)]);
   if (icnt == 0  ||  dyp <= dyptol), break, end
end  % for

% No error encountered.

ier = iter;
return;

end  % tspsp
