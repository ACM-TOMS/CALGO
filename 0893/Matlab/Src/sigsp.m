function [sigma,dsmax,ier] = sigsp(nd,t,x,y,z,xp,yp,zp,tol,sigma)
% sigsp:  Minimum tension for convexity in a parametric curve
%
% USAGE:  [sigma,dsmax,ier] = sigsp(nd,t,x,y,z,xp,yp,zp,tol,sigma);
%
%   Given an ordered sequence of points C(I) = (X(I),Y(I)) 
% or C(I) = (X(I),Y(I),Z(I)) with associated derivative
% vectors CP(I) = (XP(I),YP(I)) or CP(I) = (XP(I),YP(I),
% ZP(I)) and knots t(I), this function determines the small-
% est (nonnegative) tension factors SIGMA such that the 
% Hermite interpolatory tension spline curve c(t) preserves 
% local convexity of the data.
%
%   For knot subinterval [t1,t2], denote the endpoint values
% and derivative vectors by c1,c2 and d1,d2, respectively,
% and let d = (c2-c1)/h for h = t2-t1.  Also, define 
% e1 = d1 X d and e2 = d X d2.  (In the case of a planar 
% curve, the cross products are z-components and e1,e2 are 
% scalars.)  The data is locally convex if <e1,e2> > 0.  In
% this case the tension factor is made sufficiently large 
% that the curve segment (when projected onto the plane with
% normal e1 or e2) has no inflection point.  The curvature
% vector is proportional to k(t) = c'(t) X c''(t), and the 
% requirement for convexity is <k(t),e1> > 0 and <k(t),e2>
% > 0 for all t in [t1,t2].  It can be shown that a sufficient
% condition is positivity at the two endpoints, and this is
% equivalent to g(sigma) > 0 for the strictly increasing
% function
%
%   g(sigma) = sigma*coshm(sigma)/sinhm(sigma) - m,
%
% where m = max{<d1 X d2,e1>/|e1|^2, <d1 X d2,e2>/<e1,e2>,
% <d1 X d2,e1>/<e1,e2>, <d1 X d2,e2>/|e2|^2, 0}.  In the
% case of a space curve, even with <e1,e2> > 0, it is 
% possible that <d1 X d2,e1> and <d1 X d2,e2> could be nega-
% tive, resulting in m = 0 and g(sigma) > 0 for sigma = 0.
% (g(0) = 3-m).
%
%   Note that the data defines a sign for discrete torsion
% on each interval:  <d1,e2> = <e1,d2> = det(d1,d,d2).  The
% interpolatory tension spline automatically preserves this
% property for all values of sigma.
%
%   SIGSP may be used in conjunction with Function YPC2
% (or YPC2P) in order to produce a C-2 interpolant that
% preserves the shape properties of the data.  This is
% achieved by calling YPC2 with SIGMA initialized to the
% zero vector, and then alternating calls to SIGSP with
% calls to YPC2 until the change in SIGMA is small (refer to
% the parameter descriptions for SIGMA, DSMAX and IER), or
% the maximum relative change in YP is bounded by a toler-
% ance (a reasonable value is .01).
%
%   Refer to Function SIGBP for a means of selecting mini-
% mum tension factors to satisfy bounds constraints.
%
% On input:
%
%       ND = Number of dimensions:
%            ND = 2 if the points lie in a plane.
%            ND = 3 if the points are in 3-space.
%
%       T = Vector of length N containing a strictly in-
%           creasing sequence of knots (discrete parameter
%           values).  Refer to Function ARCL2D or ARCL3D.
%           N >= 2.
%
%       X,Y,Z = Vectors of length N containing the Cartesian
%               coordinates of an ordered sequence of data
%               points C(I), I = 1 to N.  The curve is con-
%               strained to pass through these points:  
%               c(t(I)) = C(I).  Z is an unused dummy par-
%               ameter if ND = 2.
%
%       XP,YP,ZP = Vectors of length N containing the compo-
%                  nents of first derivative vectors CP(I)
%                  for I = 1 to N.  Refer to Functions 
%                  YPC1T, YPC2, and YPC2P.
%
%       TOL = Nonnegative tolerance for the zero finder when
%             nonzero finite tension is necessary and
%             sufficient to satisfy the constraint.  Use
%             TOL = 0 for full accuracy.
%
%       SIGMA = Vector of length N-1 containing minimum 
%               values of the tension factors.  SIGMA(I) is
%               associated with interval (t(I),t(I+1)) and 
%               SIGMA(I) >= 0 for I = 1 to N-1.  SIGMA 
%               should be set to the zero vector if minimal 
%               tension is desired, and should be unchanged
%               from a previous call in order to ensure con-
%               vergence of the C-2 iterative procedure.
%
% On output:
%
%       SIGMA = Array containing tension factors for which
%               Hermite interpolatory tension spline c
%               preserves local convexity of the data,
%               with the restriction that SIGMA(I) <= SBIG
%               for all I (unless the input value is larger).
%               The factors are as small as possible (within
%               the tolerance), but not less than their
%               input values.  In interval [t(I),t(I+1)], if
%               the data is not locally convex (or the end-
%               point values are not distinct, or either
%               endpoint derivative vector is zero), then
%               SIGMA(I) = 0 (unless the input value is pos-
%               itive), and thus c is cubic in the interval.
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
%             IER = -1 if ND or N is outside its valid
%                      range.  SIGMA is not altered in
%                      this case.
%             IER = -I if t(I) <= t(I-1) for some I in the
%                      range 2 to N.  SIGMA(J-1) is unal-
%                      tered for J = I to N in this case.
%
% Module required by SIGSP:  SNHCSH
%
%***********************************************************

global SBIG

% Set fid = 1 to print diagnostic error messages.

fid = -1;

% Initialize change counter IER and maximum change DSMAX, and
%   test for ND invalid or N < 2.

ier = 0;
dsmax = 0;
n = length(t);
if (nd < 2  ||  nd > 3  ||  n < 2) 
   ier = -1;
   return;
end

% Compute an absolute tolerance FTOL = abs(TOL) and a
%   relative tolerance RTOL = 1000*MACHEPS.

ftol = abs(tol);
rtol = 1000.0*eps;

% Loop on subintervals.

for i = 1:n-1
   if (fid > 0) 
      fprintf(fid,'\n\n SIGSP:  Interval %0.0f\n', i);
   end
   ip1 = i + 1;
   dt = t(ip1) - t(i);
   if (dt <= 0)
      ier = -ip1;
      return;
   end
   sigin = sigma(i);
   if (sigin >= SBIG), continue, end

% Compute parameters for interval [t(i),t(ip1)].

   if (nd == 2)
      d = [x(ip1)-x(i) y(ip1)-y(i)]./dt;
      d1 = [xp(i) yp(i)];
      d2 = [xp(ip1) yp(ip1)];
      e1 = d1(1)*d(2) - d1(2)*d(1);
      e2 = d(1)*d2(2) - d(2)*d2(1);
      d1cd2 = d1(1)*d2(2) - d1(2)*d2(1);
   else
      d = [x(ip1)-x(i) y(ip1)-y(i) z(ip1)-z(i)]./dt;
      d1 = [xp(i) yp(i) zp(i)];
      d2 = [xp(ip1) yp(ip1) zp(ip1)];
      e1 = cross(d1,d);
      e2 = cross(d,d2);
      d1cd2 = cross(d1,d2);
   end
   e1e2 = dot(e1,e2);
   e1e1 = dot(e1,e1);
   e2e2 = dot(e2,e2);
   dde1 = dot(d1cd2,e1);
   dde2 = dot(d1cd2,e2);

   while (true)

% Test for SIGMA = 0 sufficient.

      sig = 0;
      if (e1e2 <= 0), break, end

      m = max([dde1/e1e1,dde2/e1e2,dde1/e1e2,dde2/e2e2,0]);
      if (3.0 - m >= 0), break, end

% Convexity:  Find a zero of F(SIG) = SIG*COSHM(SIG)/
%   SINHM(SIG) - m.
%
%   F(0) = 3-m < 0, F(m) >= 0, the derivative of F
%     vanishes at SIG = 0, and the second derivative of F is
%     .2 at SIG = 0.  A quadratic approximation is used to
%     obtain a starting point for the Newton method.

      sig = sqrt(10.0*m-30.0);
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

         f = sig*t1 - m;
         if (fid > 0) 
            fprintf(fid,['                 SIG = %15.8e, F(SIG) = ', ...
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

end  % sigsp
