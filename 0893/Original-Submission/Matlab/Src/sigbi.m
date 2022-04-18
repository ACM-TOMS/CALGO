function [sigma,icflg,dsmax,ier] = sigbi(x,y,yp,tol,b,bmax,sigma)
% sigbi:  Minimum tension for constrained interpolant
%
% USAGE:  [sigma,icflg,dsmax,ier] = sigbi(x,y,yp,tol,b,bmax,sigma);
%
%   Given a set of abscissae X with associated data values Y
% and derivatives YP, this function determines the small-
% est (nonnegative) tension factors SIGMA such that the Her-
% mite interpolatory tension spline H(x) satisfies a set of
% user-specified constraints.
%
%   SIGBI may be used in conjunction with Function YPC2
% (or YPC2P) in order to produce a C-2 interpolant that
% satisfies the constraints.  This is achieved by calling
% YPC2 with SIGMA initialized to the zero vector, and then
% alternating calls to SIGBI with calls to YPC2 until the
% change in SIGMA is small (refer to the parameter descrip-
% tions for SIGMA, DSMAX and IER), or the maximum relative
% change in YP is bounded by a tolerance (a reasonable value
% is .01).  A similar procedure may be used to produce a C-2
% shape-preserving smoothing curve (Function SMCRV).
%
%   Refer to Function SIGS for a means of selecting mini-
% mum tension factors to preserve shape properties of the
% data.
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
%       TOL = Tolerance whose magnitude determines how close
%             each tension factor is to its optimal value
%             when nonzero finite tension is necessary and
%             sufficient to satisfy a constraint.  Refer to
%             functions SIG0, SIG1, and SIG2.  TOL should be
%             set to 0 for optimal tension.
%
%       B = Array dimensioned 5 by N-1 containing bounds or
%           flags which define the constraints.  For I = 1
%           to N-1, column I defines the constraints associ-
%           ated with interval I (X(I),X(I+1)) as follows:
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
%              used as an upper bound in B (or when its
%              negative is used as a lower bound), specifies
%              that no constraint is to be enforced.
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
%               H(x) satisfies the constraints defined by B,
%               with the restriction that SIGMA(I) <= SBIG
%               for all I (unless the input value is larger).
%               The factors are as small as possible (within
%               the tolerance), but not less than their
%               input values.  If infinite tension is re-
%               quired in interval (X(I),X(I+1)), then
%               SIGMA(I) = SBIG (and H is an approximation 
%               to the linear interpolant on the interval),
%               and if no constraint is specified in the
%               interval, then SIGMA(I) = 0 (unless the
%               input value is positive), and thus H is
%               cubic.  Invalid constraints are treated as
%               null constraints.
%
%       ICFLG = Array of size 1 by N-1 containing invalid 
%               constraint flags associated with intervals.
%               For I = 1 to N-1, ICFLG(I) is a 5-bit value
%               b5b4b3b2b1, where bK = 1 if and only if 
%               constraint K cannot be satisfied.  Thus, all
%               constraints in interval I are satisfied if 
%               and only if ICFLG(I) = 0 (and IER >= 0).
%
%       DSMAX = Maximum increase in a component of SIGMA
%               from its input value.  The increase is a
%               relative change if the input value is
%               positive, and an absolute change otherwise.
%
%       IER = Error indicator and information flag:
%             IER = I if no errors (other than invalid con-
%                     straints) were encountered and I
%                     components of SIGMA were altered from
%                     their input values for 0 <= I <=
%                     N-1.
%             IER = -1 if N < 2.  SIGMA is not altered, and 
%                      ICFLG = 0 in this case.
%             IER = -I if X(I) <= X(I-1) for some I in the
%                      range 2 to N.  SIGMA(J) is not 
%                      altered, and ICFLG(J) = 0 for J >=
%                      I-1 in this case.
%
% Modules required by SIGBI:  SIG0, SIG1, SIG2, SNHCSH
%
%***********************************************************

global SBIG

% Initialize change counter IER and maximum change DSMAX, and
%   test for n < 2.

ier = 0;
dsmax = 0;
n = length(x);
icflg = zeros(1,n-1);
if (n < 2) 
   ier = -1;
   return;
end

% Loop on subintervals.

for i = 1:n-1
   if (x(i) >= x(i+1))
      ier = -(i+1);
      return;
   end

% Loop on constraints for interval I.  SIG is set to the
%   largest tension factor required to satisfy all five
%   constraints.  ICFK = 2**(K-1) is the increment for
%   ICFLG(I) when constraint K is invalid.

   sig = 0;
   icfk = 0.5;
   for k = 1:5
      icfk = 2*icfk;
      bnd = b(k,i);
      if (k < 5  &&  abs(bnd) >= bmax), continue, end
      if (k <= 2)
         ifl = 3 - 2*k;
         [s,ierr] = sig0(x(i),x(i+1),y(i),y(i+1),yp(i),yp(i+1), ...
                         ifl,bnd,tol);
      elseif (k <= 4)
         ifl = 7 - 2*k;
         [s,ierr] = sig1(x(i),x(i+1),y(i),y(i+1),yp(i),yp(i+1), ...
                         ifl,bnd,tol);
      else
         if (bnd == 0), continue, end
         ifl = -1;
         if (bnd > 0), ifl = 1; end
         [s,ierr] = sig2(x(i),x(i+1),y(i),y(i+1),yp(i),yp(i+1), ...
                         ifl,tol);
      end
      if (ierr == -2) 

%   An invalid constraint was encountered.  Increment ICFLG(I).

         icflg(i) = icflg(i) + icfk;
      else

%   Update SIG.

         sig = max([sig,s]);
      end

%   Bottom of loop on constraints K.

   end

% Bottom of loop on intervals:  update SIGMA(I), IER, and
%   DSMAX if necessary.

   sig = min([sig,SBIG]);
   sigin = sigma(i);
   if (sig > sigin) 
      sigma(i) = sig;
      ier = ier + 1;
      dsig = sig-sigin;
      if (sigin > 0), dsig = dsig/sigin; end
      dsmax = max([dsmax,dsig]);
   end
end

% No errors (other than invalid constraints) encountered.

return;

end  % sigbi
