function [sig,ier] = sig2(x1,x2,y1,y2,y1p,y2p,ifl,tol)
% sig2:  Minimum tension factor for convexity
%
% USAGE:  [sig,ier] = sig2(x1,x2,y1,y2,y1p,y2p,ifl,tol);
%
%   Given a pair of abscissae with associated ordinates and
% slopes, this function determines the smallest (nonnega-
% tive) tension factor SIG such that the Hermite interpo-
% latory tension spline H(x) preserves convexity (or con-
% cavity) of the data;  i.e.,
%
%   Y1P <= S <= Y2P implies HPP(x) >= 0  or
%   Y1P >= S >= Y2P implies HPP(x) <= 0
%
% for all x in the open interval (X1,X2), where S = (Y2-Y1)/
% (X2-X1) and HPP denotes the second derivative of H.  Note,
% however, that infinite tension is required if Y1P = S or
% Y2P = S (unless Y1P = Y2P = S).
%
% On input:
%
%       X1,X2 = Abscissae.  X1 < X2.
%
%       Y1,Y2 = Values of H at X1 and X2.
%
%       Y1P,Y2P = Derivative values of H at X1 and X2.
%
%       IFL = Option indicator (sign of HPP):
%             IFL = -1 if HPP is to be bounded above by 0.
%             IFL = 1 if HPP is to be bounded below by 0
%                     (preserve convexity of the data).
%
%       TOL = Tolerance whose magnitude determines how close
%             SIG is to its optimal value when nonzero
%             finite tension is necessary and sufficient to
%             satisfy convexity or concavity.  In the case
%             of convexity, SIG is chosen so that 0 <=
%             HPPMIN <= abs(TOL), where HPPMIN is the min-
%             imum value of HPP in the interval.  In the
%             case of concavity, the maximum value of HPP
%             satisfies -abs(TOL) <= HPPMAX <= 0.  Thus,
%             the constraint is satisfied but possibly with
%             more tension than necessary.
%
% On output:
%
%       SIG = Tension factor defined above unless IER < 0,
%             in which case SIG = -1.  If IER = 1, SIG
%             is set to SBIG, resulting in an approximation
%             to the linear interpolant of the endpoint
%             values.  Note, however, that SIG may be
%             larger than SBIG if IER = 0.
%
%       IER = Error indicator:
%             IER = 0 if no errors were encountered and fin-
%                     ite tension is sufficient to satisfy
%                     the constraint.
%             IER = 1 if no errors were encountered but in-
%                     finite tension is required to satisfy
%                     the constraint.
%             IER = -1 if X2 <= X1 or abs(IFL) ~= 1.
%             IER = -2 if the constraint cannot be satis-
%                      fied:  the sign of S-Y1P or Y2P-S
%                      does not agree with IFL.
%
% Module required by SIG2:  SNHCSH
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

% Compute the slope and second differences, and test for
%   an invalid constraint.

s = (y2-y1)/dx;
d1 = s - y1p;
d2 = y2p - s;
if ((ifl > 0  &&  min([d1,d2]) < 0)  || ...
    (ifl < 0  &&  max([d1,d2]) > 0))
   sig = -1.0;
   ier = -2;
   return;
end

% Test for infinite tension required.

if (d1*d2 == 0  &&  d1 ~= d2)
   sig = SBIG;
   ier = 1;
   return;
end

% Test for SIG = 0 sufficient.

sig = 0;
ier = 0;
if (d1*d2 == 0), return; end
t = max([d1/d2,d2/d1]);
if (t <= 2.0), return; end

% Find a zero of F(SIG) = SIG*COSHM(SIG)/SINHM(SIG) - (T+1).
%   Since the derivative of F vanishes at the origin, a
%   quadratic approximation is used to obtain an initial
%   estimate for the Newton method.

tp1 = t + 1.0;
sig = sqrt(10.0*t-20.0);
nit = 0;
if fid > 0
   fprintf(fid,'\n\n SIG2:  F(0) = %15.8e\n', tp1);
end

%   Compute an absolute tolerance FTOL = abs(TOL) and a
%     relative tolerance RTOL = 1000*MACHEPS.

ftol = abs(tol);
rtol = 1000.0*eps;

% Top of loop.

while (true)
   if (sig <= 0.5) 

%   Evaluate F and its derivative FP.
%   Use approximations designed to avoid cancellation error
%     in the hyperbolic functions.

      [sinhm,coshm] = snhcsh(sig);
      t1 = coshm/sinhm;
      fp = t1 + sig*(sig/sinhm - t1*t1 + 1.0);
   else

%   Scale SINHM and COSHM by 2*exp(-SIG) in order to avoid
%     overflow.

      ems = exp(-sig);
      ssm = 1.0 - ems*(ems+sig+sig);
      t1 = (1.0-ems)*(1.0-ems)/ssm;
      fp = t1 + sig*(2.0*sig*ems/ssm - t1*t1 + 1.0);
   end

   f = sig*t1 - tp1;
   nit = nit + 1;
   if (fid > 0) 
      fprintf(fid,'    %0.0f:  SIG = %15.8e, F = %15.8e\n', nit, sig, f);
   end

%   Test for convergence.

   if (fp <= 0), return; end
   dsig = -f/fp;
   if (abs(dsig) <= rtol*sig  ||  (f >= 0  &&  f <= ftol)  ||  ...
       abs(f) <= rtol), return; end

%   Update SIG.

   sig = sig + dsig;
end  % while

end  % sig2
