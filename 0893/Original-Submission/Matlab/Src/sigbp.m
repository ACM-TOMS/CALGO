function [sigma,dsmax,ier] = sigbp(x,y,xp,yp,tol,bl,bu,bmax,sigma)
% sigbp:  Minimum tension for constrained planar curve
%
% USAGE:  [sigma,dsmax,ier] = sigbp(x,y,xp,yp,tol,bl,bu,bmax,sigma);
%
%   Given an ordered sequence of points C(I) = (X(I),Y(I))
% with associated derivative vectors CP(I) = (XP(I),YP(I)),
% this function determines the smallest (nonnegative) ten-
% sion factors SIGMA such that a parametric planar curve
% C(t) satisfies a set of user-specified constraints.  The
% components x(t) and y(t) of C(t) are the Hermite interpo-
% latory tension splines defined by the data and tension
% factors:  C(t(I)) = C(I) and C'(t(I)) = CP(I) for para-
% meter values t(1), t(2), ..., t(N).  In each subinterval
% [t1,t2], the signed perpendicular distance from the
% corresponding line segment C1-C2 to the curve C(t) is
% given by the vector cross product
%
%     d(t) = (C2-C1)/DC X (C(t)-C1)
%
% where DC = abs(C2-C1) is the length of the line segment.
% The associated tension factor SIGMA is chosen to satisfy
% an upper bound on the maximum of d(t) and a lower bound on
% the minimum of d(t) over t in [t1,t2].  Thus, the upper
% bound is associated with distance to the left of the line
% segment as viewed from C1 toward C2.  Note that the curve
% is assumed to be parameterized by arc length (Function
% ARCL2D) so that t2-t1 = DC.  If this is not the case, the
% required bounds should be scaled by DC/(t2-t1) to obtain
% the input parameters BL and BU.
%
%   SIGBP may be used in conjunction with Function YPC2
% (or YPC2P) in order to produce a C-2 interpolant that
% satisfies the constraints.  This is achieved by calling
% YPC2 with SIGMA initialized to the zero vector, and then
% alternating calls to SIGBP with calls to YPC2 until the
% change in SIGMA is small (refer to the parameter descrip-
% tions for SIGMA, DSMAX and IER), or the maximum relative
% change in YP is bounded by a tolerance (a reasonable value
% is .01).
%
% On input:
%
%       X,Y = Vectors of length N containing the Cartesian
%             coordinates of the points C(I), I = 1 to N.
%             N >= 2.
%
%       XP,YP = Vectors of length N containing the components
%               of the derivative (velocity) vectors CP(I).
%               Refer to Functions YPC1, YPC1P, YPC2,
%               YPC2P, and SMCRV.
%
%       TOL = Nonnegative tolerance for the zero finder when
%             nonzero finite tension is necessary and 
%             sufficient to satisfy the constraint.  Use
%             TOL = 0 for full accuracy.
%
%       BL,BU = Vectors of length N-1 containing lower and
%               upper bounds, respectively, which define
%               the constraints as described above.  BL(I)
%               < 0 and BU(I) > 0 for I = 1 to N-1.  A null
%               straint is specified by BL(I) <= -BMAX or
%               BU(I) >= BMAX.
%
%       BMAX = User-defined value of infinity which, when
%              used as an upper bound in BU (or when its
%              negative is used as a lower bound in BL),
%              specifies that no constraint is to be en-
%              forced.
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
%               d(t) satisfies the constraints defined by
%               BL and BU, with the restriction that
%               SIGMA(I) <= SBIG for all I (unless the input
%               value is larger).  The factors are as small
%               as possible (within the tolerance), but not
%               less than their input values.  If no con-
%               straint is specified in interval I, then
%               SIGMA(I) = 0 (unless the input value is
%               positive), and thus x(t) and y(t) are cubic
%               polynomials.
%
%       DSMAX = Maximum increase in a component of SIGMA
%               from its input value.  The increase is a
%               relative change if the input value is
%               positive, and an absolute change otherwise.
%
%       IER = Error indicator and information flag:
%             IER = I if no errors were encountered and I
%                     components of SIGMA were altered from
%                     their input values for 0 <= I <=
%                     N-1.
%             IER = -1 if N < 2.  SIGMA is not altered in
%                      this case.
%             IER = -I if BL(I-1) >= 0 or BU(I-1) <= 0
%                      for some I in the range 2 to N.
%                      SIGMA(J) is unaltered for J >= I-1
%                      in this case.
%
% Module required by SIGBP:  SNHCSH
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

% Compute options for fzero.

tol = max(tol,eps);
options = optimset('TolX',tol);

% Loop on subintervals.

for i = 1:n-1
   ip1 = i + 1;
   blo = bl(i);
   bhi = bu(i);
   sigin = sigma(i);
   if (fid > 0) 
      fprintf(fid,['\n\n SIGBP:  Interval %0.0f, BL = %10.3e, ', ...
              'BU = %10.3e, SIGIN = %15.8e\n'], i, blo, bhi, sigin);
   end
   if (blo >= 0  ||  bhi <= 0) 
      ier = -(ip1);
      return;
   end
   if (sigin >= SBIG), continue, end

% Initialize SIG to 0 and test for a null constraint.

   sig = 0;
   if (blo <= -bmax  &&  bhi >= bmax) 

% Update SIGMA(I), IER, and DSMAX if necessary, and continue
%   to the next subinterval.

      if (sig > sigin) 
         sigma(i) = sig;
         ier = ier + 1;
         dsig = sig-sigin;
         if (sigin > 0), dsig = dsig/sigin; end
         dsmax = max([dsmax,dsig]);
      end
      continue
   end

% Test for SIG = 0 sufficient.
%
%   The signed orthogonal distance is d(b) = b*(1-b)*
%     (b*V1 - (1-b)*V2), where b = (t2-t)/(t2-t1),
%     V1 = (C2-C1) X CP(1), and V2 = (C2-C1) X CP(2).

   dx = x(ip1) - x(i);
   dy = y(ip1) - y(i);
   v1 = dx*yp(i) - dy*xp(i);
   v2 = dx*yp(ip1) - dy*xp(ip1);

%   Set DP and DM to the maximum and minimum values of d(b)
%     for b in [0,1].  Note that DP >= 0 and DM <= 0.

   s = v1 + v2;
   if (s == 0) 

%   The derivative d'(b) is zero at the midpoint b = .5.

      if (v1 >= 0) 
         dp = v1/4.0;
         dm = 0;
      else
         dp = 0;
         dm = v1/4.0;
      end
   else

%   Set RP/RM to the roots of the quadratic equation d'(b) =
%     (B0 +/- SQRT(D0))/(3*S) = V2/(B0 -/+ SQRT(D0)) = 0,
%     where B0 = V1 + 2*V2 and D0 = V1^2 + V1*V2 + V2^2.
%     The expression is chosen to avoid cancellation error.

      b0 = s + v2;
      d0 = s*s - v1*v2;
      t = b0 + sign(b0)*sqrt(d0);
      if (b0 >= 0) 
         rp = t/(3.0*s);
         rm = v2/t;
      else
         rp = v2/t;
         rm = t/(3.0*s);
      end
      if (v1 <= 0  &&  v2 >= 0) 

%   The maximum is DP = 0 at the endpoints.

         dp = 0;
      else
         dp = rp*(1.0-rp)*(rp*s - v2);
      end
      if (v1 >= 0  &&  v2 <= 0)

%   The minimum is DM = 0 at the endpoints.

         dm = 0;
      else
         dm = rm*(1.0-rm)*(rm*s - v2);
      end
   end

%   SIG = 0 is sufficient to satisfy the constraints iff
%     DP <= BHI and DM >= BLO iff F0 >= 0.

   f0 = min([bhi-dp, dm-blo]);
   if (f0 >= 0) 

% Update SIGMA(I), IER, and DSMAX if necessary, and continue
%   to the next subinterval.

      if (sig > sigin) 
         sigma(i) = sig;
         ier = ier + 1;
         dsig = sig-sigin;
         if (sigin > 0), dsig = dsig/sigin; end
         dsmax = max([dsmax,dsig]);
      end
      continue
   end

% Find a zero of F(SIG) = min(BHI-DP,DM-BLO), where DP and
%   DM are the maximum and minimum values of d(b).  F is an
%   increasing function, F(0) = F0 < 0, and F = FMAX =
%   min(BHI,-BLO) for SIG sufficiently large.  Note that F
%   has a discontinuity in its first derivative if the
%   curves BHI-DP and DM-BLO (as functions of SIG) inter-
%   sect, and the rate of convergence of the zero finder is
%   reduced to linear if such an intersection occurs near
%   the zero of F.
%
% Store shared variables needed by nested function fsigbp.

   fmax = min([bhi,-blo]);
   v2m1 = v2 - v1;
   nit = -1;

   f = fsigbp(SBIG);

   if (fid > 0) 
      fprintf(fid,[' F(0) = %15.8e, F(SBIG) = %15.8e, ', ...
                   'FMAX = %15.8e\n\n'], f0, f, fmax);
   end
   if f <= 0
      sig = SBIG;
   else

% [0,SBIG] is a bracketing interval.

      nit = 0;
      sig = fzero(@fsigbp,[0,SBIG],options);
   end

% Bottom of loop on intervals:  update SIGMA(I), IER, and
%   DSMAX if necessary.

   if (sig > sigin) 
      sigma(i) = sig;
      ier = ier + 1;
      dsig = sig-sigin;
      if (sigin > 0), dsig = dsig/sigin; end
      dsmax = max([dsmax,dsig]);
   end
end

% No errors encountered.

return;

   function f = fsigbp(sig)
% Nested function for evaluation of F.

      if sig == 0
         f = f0;
         return;
      end
      ems = exp(-sig);
      if (sig <= 0.5)

%   SIG <= .5:  use approximations designed to avoid can-
%                 cellation error (associated with small
%                 SIG) in the modified hyperbolic functions.

         [sinhm,coshm,coshmm] = snhcsh(sig);
         sinh = sinhm + sig;
         a1 = sig*coshm*v2 - sinhm*v2m1;
         a2 = sig*sinh*v2 - coshm*v2m1;
         a = a2 - a1;
         aa = a/ems;
         e = sig*sinhm - coshmm - coshmm;
      else

%   SIG > .5:  scale SINHM and COSHM by 2*EXP(-SIG) in order
%              to avoid overflow.

         tm = 1.0 - ems;
         sinh = tm*(1.0+ems);
         sinhm = sinh - 2.0*sig*ems;
         coshm = tm*tm;
         a1 = sig*coshm*v2 - sinhm*v2m1;
         a2 = sig*sinh*v2 - coshm*v2m1;
         aa = 2.0*(sig*tm*v2 + (tm-sig)*v2m1);
         a = ems*aa;
         e = sig*sinh - coshm - coshm;
      end
      if (s == 0) 

%   The derivative d'(b) is zero at the midpoint b = .5.

         eb = sig*coshm - sinhm - sinhm;
         if (v1 >= 0) 
            dp = e*v1/(sig*(sqrt(eb*eb-e*e)+eb));
            dm = 0;
         else
            dp = 0;
            dm = e*v1/(sig*(sqrt(eb*eb-e*e)+eb));
         end
         f = min([bhi-dp, dm-blo]);
      else

%   d'(b)*DC = V2 - (A1*sinh(SIG*b) - A2*coshm(SIG*b))/E = 0
%     for ESB = (-B +/- sqrt(D))/A = C/(-B -/+ sqrt(D)),
%     where ESB = exp(SIG*b), A = A2-A1, D = B^2 - A*C, and
%     B and C are defined below.

         b = -coshm*s;
         c = a2 + a1;
         d = b*b - a*c;
         if (d < 0)
            f = fmax;
         else
            t1 = sqrt(d);
            t = -b - sign(b)*t1;

            rsp = 0;
            if (b < 0  &&  aa ~= 0) 
               if (t/aa > 0), rsp = sig + log(t/aa); end
            end
            if ((b > 0  ||  aa == 0)  &&  c/t > 0), rsp = log(c/t); end
            if ((rsp <= 0  ||  rsp >= sig)  &&  b ~= 0)

%   The maximum is DP = 0 at the endpoints.

               dp = 0;
            else
               dp = -(b*rsp+a1+t1)/(sig*e);
            end

            rsm = 0;
            if (b > 0  &&  aa ~= 0)
               if (t/aa > 0), rsm = sig + log(t/aa); end
            end
            if ((b < 0  ||  aa == 0)  &&  c/t > 0), rsm = log(c/t); end
            if ((rsm <= 0  ||  rsm >= sig)  &&  b ~= 0) 

%   The minimum is DM = 0 at the endpoints.

               dm = 0;
            else
               dm = -(b*rsm+a1-t1)/(sig*e);
            end
            f = min([bhi-dp, dm-blo]);
         end
      end

%   Update the number of iterations NIT.

      nit = nit + 1;
      if (fid > 0  &&  nit > 0) 
         fprintf(fid,'    %0.0f:  SIG = %15.8e, F = %15.8e\n', nit, sig, f);
      end
      return;
      end

end  % sigbp
