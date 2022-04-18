function [sigma,dsmax,ier] = sigs(x,y,yp,tol,sigma)
% sigs:  Minimum tension for monotonicity and convexity
%
% USAGE:  [sigma,dsmax,ier] = sigs(x,y,yp,tol,sigma);
%
%   Given a set of abscissae X with associated data values Y
% and derivatives YP, this function determines the small-
% est (nonnegative) tension factors SIGMA such that the Her-
% mite interpolatory tension spline H(x) preserves local
% shape properties of the data.  In an interval (X1,X2) with
% data values Y1,Y2 and derivatives YP1,YP2, the properties
% of the data are
%
%       Monotonicity:  S, YP1, and YP2 are nonnegative or
%                        nonpositive,
%  and
%       Convexity:     YP1 <= S <= YP2  or  YP1 >= S
%                        >= YP2,
%
% where S = (Y2-Y1)/(X2-X1).  The corresponding properties
% of H are constant sign of the first and second deriva-
% tives, respectively.  Note that, unless YP1 = S = YP2, in-
% finite tension is required (and H is linear on the inter-
% val) if S = 0 in the case of monotonicity, or if YP1 = S
% or YP2 = S in the case of convexity.
%
%   SIGS may be used in conjunction with Function YPC2
% (or YPC2P) in order to produce a C-2 interpolant that
% preserves the shape properties of the data.  This is
% achieved by calling YPC2 with SIGMA initialized to the
% zero vector, and then alternating calls to SIGS with
% calls to YPC2 until the change in SIGMA is small (refer to
% the parameter descriptions for SIGMA, DSMAX and IER), or
% the maximum relative change in YP is bounded by a toler-
% ance (a reasonable value is .01).  A similar procedure may
% be used to produce a C-2 shape-preserving smoothing curve
% (Function SMCRV).
%
%   Refer to Function SIGBI for a means of selecting mini-
% mum tension factors to satisfy more general constraints.
%
% On input:
%
%       X = Vector of length N containing a strictly in-
%           creasing sequence of abscissae:  X(I) < X(I+1)
%           for I = 1 to N-1.  N >= 2.
%
%       Y = Vector of length N containing data values (or
%           function values computed by SMCRV) associated
%           with the abscissae.  H(X(I)) = Y(I) for I =
%           1 to N.
%
%       YP = Vector of length N containing first derivatives
%            of H at the abscissae.  Refer to Functions
%            YPC1, YPC1P, YPC2, YPC2P, and SMCRV.
%
%       TOL = Nonnegative tolerance for the zero finder when
%             nonzero finite tension is necessary and
%             sufficient to satisfy the constraint.  Use
%             TOL = 0 for full accuracy.
%
%       SIGMA = Vector of length N-1 containing minimum val-
%               ues of the tension factors.  SIGMA(I) is as-
%               sociated with interval (I,I+1) and SIGMA(I)
%               >= 0 for I = 1 to N-1.  SIGMA should be
%               set to the zero vector if minimal tension
%               is desired, and should be unchanged from a
%               previous call in order to ensure convergence
%               of the C-2 iterative procedure.
%
% On output:
%
%       SIGMA = Array containing tension factors for which
%               H(x) preserves the properties of the data,
%               with the restriction that SIGMA(I) <= SBIG
%               for all I (unless the input value is larger).
%               The factors are as small as possible (within
%               the tolerance), but not less than their
%               input values.  If infinite tension is re-
%               quired in interval (X(I),X(I+1)), then
%               SIGMA(I) = SBIG (and H is an approximation
%               to the linear interpolant on the interval),
%               and if neither property is satisfied by the
%               data, then SIGMA(I) = 0 (unless the input
%               value is positive), and thus H is cubic in
%               the interval.
%
%       DSMAX = Maximum increase in a component of SIGMA
%               from its input value.  The increase is a
%               relative change if the input value is
%               nonzero, and an absolute change otherwise.
%
%       IER = Error indicator and information flag:
%             IER = I if no errors were encountered and I
%                     components of SIGMA were altered from
%                     their input values for 0 <= I <=
%                     N-1.
%             IER = -1 if N < 2.  SIGMA is not altered in
%                      this case.
%             IER = -I if X(I) <= X(I-1) for some I in the
%                      range 2 to N.  SIGMA(J-1) is unal-
%                      tered for J = I to N in this case.
%
% Module required by SIGS:  SNHCSH
%
%***********************************************************

global SBIG

% Set fid = 1 to print diagnostic error messages.

fid = -1;

% Initialize change counter IER and maximum change DSMAX, and
%   test for n < 2.

ier = 0;
dsmax = 0;
n = length(x);
if (n < 2) 
   ier = -1;
   return;
end

% Compute an absolute tolerance FTOL = abs(TOL) and a
%   relative tolerance RTOL = 1000*MACHEPS for the 
%   convexity computation, and compute options for fzero
%   (used for the monotonicity computation).

ftol = abs(tol);
rtol = 1000.0*eps;
options = optimset('TolX',max(tol,eps));

% Loop on subintervals.

for i = 1:n-1
   if (fid > 0) 
      fprintf(fid,'\n\n SIGS:  Interval %0.0f\n', i);
   end
   ip1 = i + 1;
   dx = x(ip1) - x(i);
   if (dx <= 0)
      ier = -ip1;
      return;
   end
   sigin = sigma(i);
   if (sigin >= SBIG), continue, end

% Compute first and second differences.

   s1 = yp(i);
   s2 = yp(ip1);
   s = (y(ip1)-y(i))/dx;
   d1 = s - s1;
   d2 = s2 - s;
   d1d2 = d1*d2;

   while (true)

% Test for infinite tension required to satisfy either
%   property.

      sig = SBIG;
      if ((d1d2 == 0  &&  s1 ~= s2)  || ...
          (s == 0  &&  s1*s2 > 0)), break, end

% Test for SIGMA = 0 sufficient.  The data satisfies convex-
%   ity iff D1D2 >= 0, and D1D2 = 0 implies S1 = S = S2.

      sig = 0;
      if (d1d2 >= 0)
         if (d1d2 == 0), break, end
         t = max([d1/d2,d2/d1]);
         if (t <= 2.0), break, end
         tp1 = t + 1.0;

% Convexity:  Find a zero of F(SIG) = SIG*COSHM(SIG)/
%   SINHM(SIG) - TP1.
%
%   F(0) = 2-T < 0, F(TP1) >= 0, the derivative of F
%     vanishes at SIG = 0, and the second derivative of F is
%     .2 at SIG = 0.  A quadratic approximation is used to
%     obtain a starting point for the Newton method.

         sig = sqrt(10.0*t-20.0);
         nit = 0;

%   Top of loop:

         while (true)
            if (sig <= 0.5) 
               [sinhm,coshm] = snhcsh(sig);
               t1 = coshm/sinhm;
               fp = t1 + sig*(sig/sinhm - t1*t1 + 1.0);
            else

%   Scale SINHM and COSHM by 2*EXP(-SIG) in order to avoid
%     overflow with large SIG.

               ems = exp(-sig);
               ssm = 1.0 - ems*(ems+sig+sig);
               t1 = (1.0-ems)*(1.0-ems)/ssm;
               fp = t1 + sig*(2.0*sig*ems/ssm - t1*t1 + 1.0);
            end

            f = sig*t1 - tp1;
            if (fid > 0) 
               fprintf(fid,['     CONVEXITY:  SIG = %15.8e, F(SIG) = ', ...
                       '%15.8e\n',repmat(' ',1,35),'FP(SIG) = %15.8e\n'], ...
                       sig, f, fp);
            end
            nit = nit + 1;

%   Test for convergence.
 
            if (fp <= 0), break, end
            dsig = -f/fp;
            if (abs(dsig) <= rtol*sig  ||  (f >= 0  &&  f <= ftol)  || ... 
                abs(f) <= rtol), break, end

%   Update SIG.

            sig = sig + dsig;
         end
         break
      end

% Convexity cannot be satisfied.  Monotonicity can be satis-
%   fied iff S1*S >= 0 and S2*S >= 0 since S ~= 0.

      if (s1*s < 0  ||  s2*s < 0), break, end
      t0 = 3.0*s - s1 - s2;
      d0 = t0*t0 - s1*s2;

% SIGMA = 0 is sufficient for monotonicity iff S*T0 >= 0
%   or D0 <= 0.

      if (d0 <= 0  ||  s*t0 >= 0), break, end

% Monotonicity:  find a zero of F(SIG) = SIGN(S)*HP(R),
%   where HPP(R) = 0 and HP, HPP denote derivatives of H.
%   F has a unique zero, F(0) = F0 < 0, and F approaches 
%   abs(S) as SIG increases.
%
% Store shared variables needed by nested function fsigs.

      sgn = sign(s);
      f0 = sgn*d0/(2*t0-s1-s2);
      sig = SBIG;
      fmax = sgn*(sig*s-s1-s2)/(sig-2.0);
      if (fmax <= 0), break, end
      d1pd2 = d1 + d2;
      nit = -1;

      f = fsigs(SBIG);
      if (fid > 0)
         fprintf(fid,['\n     MONOTONICITY:  F(0) = %15.8e, ', ...
                      'F(SBIG) = %15.8e\n'], f0, f); 
      end
      if (f <= 0), break, end

% [0,SBIG] is a bracketing interval.

      nit = 0;
      sig = fzero(@fsigs,[0,SBIG],options);
      break;
   end % while

%  Update SIGMA(I), IER, and DSMAX if necessary.

   sig = min([sig, SBIG]);
   if (sig > sigin) 
      sigma(i) = sig;
      ier = ier + 1;
      dsig = sig-sigin;
      if (sigin > 0), dsig = dsig/sigin; end
      dsmax = max([dsmax,dsig]);
   end
end  % for

% No errors encountered.

return;

         function f = fsigs(sig)
% Nested function for evaluation of F.

         if sig == 0
            f = f0;
            return;
         end
         if (sig <= 0.5) 

%   Use approximations to the hyperbolic functions designed
%     to avoid cancellation error with small SIG.

            [sinhm,coshm,coshmm] = snhcsh(sig);
            c1 = sig*coshm*d2 - sinhm*d1pd2;
            c2 = sig*(sinhm+sig)*d2 - coshm*d1pd2;
            a = c2 - c1;
            e = sig*sinhm - coshmm - coshmm;
            f = (sgn*(e*s2-c2) + sqrt(a*(c2+c1)))/e;
         else

%   Scale SINHM and COSHM by 2*EXP(-SIG) in order to avoid
%     overflow with large SIG.

            ems = exp(-sig);
            ems2 = ems + ems;
            tm = 1.0 - ems;
            ssinh = tm*(1.0+ems);
            ssm = ssinh - sig*ems2;
            scm = tm*tm;
            c1 = sig*scm*d2 - ssm*d1pd2;
            c2 = sig*ssinh*d2 - scm*d1pd2;

%   R is in (0,1) and well-defined iff HPP(X1)*HPP(X2) < 0.

            f = fmax;
            if (c1*(sig*scm*d1 - ssm*d1pd2) < 0)
               a = ems2*(sig*tm*d2 + (tm-sig)*d1pd2);
               if (a*(c2+c1) >= 0)
                  e = sig*ssinh - scm - scm;
                  f = (sgn*(e*s2-c2) + sqrt(a*(c2+c1)))/e;
               end
            end
         end

%   Update number of iterations NIT.

         nit = nit + 1;
         if (fid > 0  &&  nit > 0)
            fprintf(fid,[repmat(' ',1,11),'%0.0f:  SIG = %15.8e', ...
                    ', F = %15.8e\n'], nit, sig, f);
         end
         return;
         end

end  % sigs
