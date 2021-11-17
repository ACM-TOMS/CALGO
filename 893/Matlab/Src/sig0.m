function [sig,ier] = sig0(x1,x2,y1,y2,y1p,y2p,ifl,hbnd,tol)
% sig0:  Minimum tension factor for bound on function values
%
% USAGE:  [sig,ier] = sig0(x1,x2,y1,y2,y1p,y2p,ifl,hbnd,tol);
%
%   Given a pair of abscissae with associated ordinates and
% slopes, this function determines the smallest (nonnega-
% tive) tension factor SIG such that the Hermite interpo-
% latory tension spline H(x), defined by SIG and the data,
% is bounded (either above or below) by HBND for all x in
% (X1,X2).
%
% On input:
%
%       X1,X2 = Abscissae.  X1 < X2.
%
%       Y1,Y2 = Values of H at X1 and X2.
%
%       Y1P,Y2P = Derivative values of H at X1 and X2.
%
%       IFL = Option indicator:
%             IFL = -1 if HBND is a lower bound on H.
%             IFL = 1 if HBND is an upper bound on H.
%
%       HBND = Bound on H.  If IFL = -1, HBND <= min(Y1,
%              Y2).  If IFL = 1, HBND >= max(Y1,Y2).
%
%       TOL = Nonnegative tolerance for the zero finder when
%             nonzero finite tension is necessary and 
%             sufficient to satisfy the constraint.  Use
%             TOL = 0 for full accuracy.
%
% On output:
%
%       SIG = Minimum tension factor defined above unless
%             IER < 0, in which case SIG = -1.  If IER =
%             1, SIG = SBIG, resulting in an approximation
%             to the linear interpolant of the endpoint
%             values.
%
%       IER = Error indicator:
%             IER = 0 if no errors were encountered and the
%                     constraint can be satisfied with fin-
%                     ite tension.
%             IER = 1 if no errors were encountered but SIG
%                     > SBIG is required to satisfy the
%                     constraint.
%             IER = -1 if X2 <= X1 or abs(IFL) ~= 1.
%             IER = -2 if HBND is outside its valid range
%                      on input.
%
% Module required by SIG0:  SNHCSH
%
%***********************************************************

global SBIG

% Set fid = 1 to print diagnostic error messages.

fid = -1;

% Test for error 1.

dx = x2 - x1;
if (abs(ifl) ~= 1.0  ||  dx <= 0)
   sig = -1.0;
   ier = -1;
   return;
end

% Test for a valid constraint.

if ( (ifl < 0  &&  min([y1,y2]) < hbnd)  || ...
     (ifl > 0  &&  hbnd < max([y1,y2])) )
   sig = -1.0;
   ier = -2;
   return;
end

% Test for infinite tension required.

s1 = y1p;
s2 = y2p;
if ((y1 == hbnd  &&  ifl*s1 > 0)  || ...
    (y2 == hbnd  &&  ifl*s2 < 0)) 
   sig = SBIG;
   ier = 1;
   return;
end

% Test for SIG = 0 sufficient.

sig = 0;
ier = 0;
if (ifl*s1 <= 0  &&  ifl*s2 >= 0), return; end

% Compute coefficients A0 and B0 of the Hermite cubic in-
%   terpolant H0(x) = Y2 - DX*(S2*R + B0*R^2 + A0*R^3/3)
%   where R = (X2-x)/DX.

s = (y2-y1)/dx;
t0 = 3.0*s - s1 - s2;
a0 = 3.0*(s-t0);
b0 = t0 - s2;
d0 = t0*t0 - s1*s2;

% H0 has local extrema in (X1,X2) iff S1*S2 < 0 or
%   (T0*(S1+S2) < 0 and D0 >= 0).

if (s1*s2 >= 0  &&  (t0*(s1+s2) >= 0  ||  d0 < 0)), return; end
if (a0 == 0) 

% H0 is quadratic and has an extremum at R = -S2/(2*B0).
%   H0(R) = Y2 + DX*S2^2/(4*B0).  Note that A0 = 0 im-
%   plies 2*B0 = S1-S2, and S1*S2 < 0 implies B0 ~= 0.
%   Also, the extremum is a min iff HBND is a lower bound.

   f0 = (hbnd - y2 - dx*s2*s2/(4.0*b0))*ifl;
else

% A0 ~= 0 and H0 has extrema at R = (-B0 +/- SQRT(D0))/
%   A0 = S2/(-B0 -/+ SQRT(D0)), where the negative root
%   corresponds to a min.  The expression for R is chosen
%   to avoid cancellation error.  H0(R) = Y2 + DX*(S2*B0 +
%   2*D0*R)/(3*A0).

   t = -b0 - sign(b0)*sqrt(d0);
   r = t/a0;
   if (ifl*b0 > 0), r = s2/t; end
   f0 = (hbnd - y2 - dx*(s2*b0+2.0*d0*r)/(3.0*a0))*ifl;
end

%   F0 >= 0 iff SIG = 0 is sufficient to satisfy the
%     constraint.

if (f0 >= 0), return; end

% Find a zero of F = FSIG0(SIG) = (BND-H(R))*RF where the 
%   derivative of H, HP, vanishes at R.  F is generally a 
%   nondecreasing function with F(0) < 0 and F = FMAX for 
%   SIG sufficiently large.
%
% Store shared variables needed by nested function fsig0.

fmax = min([abs(y1-hbnd), abs(y2-hbnd)]);
d2 = s2 - s;
d1pd2 = s2 - s1;
nit = -1;

f = fsig0(SBIG);

if (fid > 0  &&  ifl < 0)
   fprintf(fid,['\n\n SIG0 (lower bound):  F(0) = %15.8e, ', ...
                'F(SBIG) = %15.8e\n', repmat(' ',1,46), ...
                'for SBIG = %15.8e\n'], f0, f, SBIG);
elseif (fid > 0  &&  ifl > 0)
   fprintf(fid,['\n\n SIG0 (upper bound):  F(0) = %15.8e, ', ...
                'F(SBIG) = %15.8e\n', repmat(' ',1,46), ...
                'for SBIG = %15.8e\n'], f0, f, SBIG);
end
if f <= 0
   sig = SBIG;
   ier = 1;
   return;
end

% [0,SBIG] is a bracketing interval.

nit = 0;
tol = max(tol,eps);
options = optimset('TolX',tol);
sig = fzero(@fsig0,[0 SBIG],options);
return;

   function f = fsig0(sig)
% Nested function for evaluation of F.

   if sig == 0
      f = f0;
      return;
   end
   ems = exp(-sig);
   if (sig <= 0.5) 

% SIG <= .5:  use approximations designed to avoid can-
%               cellation error (associated with small
%               SIG) in the modified hyperbolic functions.

      [sinhm,coshm,coshmm] = snhcsh(sig);
      c1 = sig*coshm*d2 - sinhm*d1pd2;
      c2 = sig*(sinhm+sig)*d2 - coshm*d1pd2;
      a = c2 - c1;
      aa = a/ems;
      e = sig*sinhm - coshmm - coshmm;
   else

% SIG > .5:  scale SINHM and COSHM by 2*EXP(-SIG) in order
%            to avoid overflow.

      tm = 1.0 - ems;
      ssinh = tm*(1.0+ems);
      ssm = ssinh - 2.0*sig*ems;
      scm = tm*tm;
      c1 = sig*scm*d2 - ssm*d1pd2;
      c2 = sig*ssinh*d2 - scm*d1pd2;
      aa = 2.0*(sig*tm*d2 + (tm-sig)*d1pd2);
      a = ems*aa;
      e = sig*ssinh - scm - scm;
   end

% HP(R) = S2 - (C1*SINH(SIG*R) - C2*COSHM(SIG*R))/E = 0
%   for ESR = (-B +/- SQRT(D))/A = C/(-B -/+ SQRT(D)),
%   where ESR = EXP(SIG*R), A = C2-C1, D = B^2 - A*C, and
%   B and C are defined below.

   b = e*s2 - c2;
   c = c2 + c1;
   d = b*b - a*c;
   f = 0;
   if (aa*c ~= 0  ||  b ~= 0), f = fmax; end
   if ((aa*c ~= 0  ||  b ~= 0)  &&  d >= 0)
      t1 = sqrt(d);
      t = -b - sign(b)*t1;
      rsig = 0;
      if (ifl*b < 0  &&  aa ~= 0) 
         if (t/aa > 0), rsig = sig + log(t/aa); end
      end
      if ((ifl*b > 0  ||  aa == 0)  &&  c/t > 0), rsig = log(c/t); end
      if ((rsig > 0  &&  rsig < sig)  ||  b == 0) 

% H(R) = Y2 - DX*(B*SIG*R + C1 + RF*SQRT(D))/(SIG*E).

         f = (hbnd - y2 + dx*(b*rsig+c1+ifl*t1)/(sig*e))*ifl;
      end
   end

% Update the number of iterations NIT.

   nit = nit + 1;
   if (fid > 0  &&  nit > 0) 
      fprintf(fid,'    %0.0f:  SIG = %15.8e, F = %15.8e\n', nit, sig, f);
   end
   return;
   end

end  % sig0
