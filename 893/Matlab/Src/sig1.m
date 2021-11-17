function [sig,ier] = sig1(x1,x2,y1,y2,y1p,y2p,ifl,hpbnd,tol)
% sig1:  Minimum tension factor for bound on first derivative
%
% USAGE:  [sig,ier] = sig1(x1,x2,y1,y2,y1p,y2p,ifl,hpbnd,tol);
%
%   Given a pair of abscissae with associated ordinates and
% slopes, this function determines the smallest (nonnega-
% tive) tension factor SIG such that the derivative HP(x)
% of the Hermite interpolatory tension spline H(x), defined
% by SIG and the data, is bounded (either above or below)
% by HPBND for all x in (X1,X2).
%
% On input:
%
%       X1,X2 = Abscissae.  X1 < X2.
%
%       Y1,Y2 = Values of H at X1 and X2.
%
%       Y1P,Y2P = Values of HP at X1 and X2.
%
%       IFL = Option indicator:
%             IFL = -1 if HPBND is a lower bound on HP.
%             IFL = 1 if HPBND is an upper bound on HP.
%
%       HPBND = Bound on HP.  If IFL = -1, HPBND <=
%               min(Y1P,Y2P,S) for S = (Y2-Y1)/(X2-X1).  If
%               IFL = 1, HPBND >= max(Y1P,Y2P,S).
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
%             IER = -2 if HPBND is outside its valid range
%                      on input.
%
% Module required by SIG1:  SNHCSH
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

s1 = y1p;
s2 = y2p;
s = (y2-y1)/dx;

% Test for a valid constraint.

if ( (ifl < 0  &&  min([s1,s2,s]) < hpbnd)  || ...
     (ifl > 0  &&  hpbnd < max([s1,s2,s])) )
   sig = -1.0;
   ier = -2;
   return;
end

% Test for infinite tension required.

if (s == hpbnd  &&  (s1 ~= s  ||  s2 ~= s))
   sig = SBIG;
   ier = 1;
   return;
end

% Test for SIG = 0 sufficient.  The Hermite cubic interpo-
%   land H0 has derivative HP0(x) = S2 + 2*B0*R + A0*R^2,
%   where R = (X2-x)/DX.

sig = 0;
ier = 0;
t0 = 3.0*s - s1 - s2;
b0 = t0 - s2;
c0 = t0 - s1;
a0 = -b0 - c0;

%   HP0(R) has an extremum (at R = -B0/A0) in (0,1) iff
%     B0*C0 > 0 and the third derivative of H0 has the
%     sign of A0.

if (b0*c0 <= 0  ||  a0*ifl > 0), return; end

%   A0*RF < 0 and HP0(R) = -D0/A0 at R = -B0/A0.

d0 = t0*t0 - s1*s2;
f0 = (hpbnd + d0/a0)*ifl;
if (f0 >= 0), return; end

% Find a zero of F = FSIG1(SIG) = (BND-HP(R))*RF, where HP 
%   has an extremum at R.  F has a unique zero, F(0) = F0
%   < 0, and F = FMAX = (BND-S)*RF > 0 for SIG sufficiently 
%   large.
%
% Store shared variables needed by nested function fsig1.

fmax = (hpbnd-s)*ifl;
d1 = s - s1;
d2 = s2 - s;
d1pd2 = d1 + d2;
nit = -1;

f = fsig1(SBIG);

if (fid > 0  &&  ifl < 0)
   fprintf(fid,['\n\n SIG1 (lower bound):  F(0) = %15.8e, ', ...
                'F(SBIG) = %15.8e\n', repmat(' ',1,46), ...
                'for SBIG = %15.8e\n'], f0, f, SBIG);
elseif (fid > 0  &&  ifl > 0)
   fprintf(fid,['\n\n SIG1 (upper bound):  F(0) = %15.8e, ', ...
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
sig = fzero(@fsig1,[0 SBIG],options);
return;

   function f = fsig1(sig)
% Nested function for evaluation of F.

   if sig == 0
      f = f0;
      return;
   end
   if (sig <= 0.5)

% Use approximations designed to avoid cancellation error
%   (associated with small SIG) in the modified hyperbolic
%   functions.

      [sinhm,coshm,coshmm] = snhcsh(sig);
      c1 = sig*coshm*d2 - sinhm*d1pd2;
      c2 = sig*(sinhm+sig)*d2 - coshm*d1pd2;
      a = c2 - c1;
      e = sig*sinhm - coshmm - coshmm;
   else

% Scale SINHM and COSHM by 2*EXP(-SIG) in order to avoid
%   overflow.

      ems = exp(-sig);
      ems2 = ems + ems;
      tm = 1.0 - ems;
      sinh = tm*(1.0+ems);
      sinhm = sinh - sig*ems2;
      coshm = tm*tm;
      c1 = sig*coshm*d2 - sinhm*d1pd2;
      c2 = sig*sinh*d2 - coshm*d1pd2;
      a = ems2*(sig*tm*d2 + (tm-sig)*d1pd2);
      e = sig*sinh - coshm - coshm;
   end

% The second derivative of H(R) has a zero at EXP(SIG*R) =
%   SQRT((C2+C1)/A) and R is in (0,1) and well-defined
%   iff HPP(X1)*HPP(X2) < 0.

   f = fmax;
   t1 = a*(c2+c1);
   if (t1 >= 0)
      if (c1*(sig*coshm*d1 - sinhm*d1pd2) < 0) 

% HP(R) = (B+SIGN(A)*SQRT(A*C))/E at the critical value
%   of R, where A = C2-C1, B = E*S2-C2, and C = C2+C1.
%   NOTE THAT RF*A < 0.

         f = (hpbnd - (e*s2-c2 - ifl*sqrt(t1))/e)*ifl;
      end
   end

% Update the number of iterations NIT.

   nit = nit + 1;
   if (fid > 0  &&  nit > 0) 
      fprintf(fid,'    %0.0f:  SIG = %15.8e, F = %15.8e\n', nit, sig, f);
   end
   return;
   end

end  % sig1
