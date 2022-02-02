C
C      DNSPLIN1:  Discrete Nonlinear Spline Interpolation
C
C                         02/18/03
C
C
C   Given a sequence of n data ordinates D1, D2, ..., Dn
C associated with strictly increasing abscissae T1, T2, ...,
C Tn, n >= 2, this program constructs a discrete approxima-
C tion to a nonlinear interpolating spline function y:  the
C graph of y is a continuous curve with continuous tangent
C and piecewise continuous curvature (with discontinuities
C allowed only at the knots (Ti,Di)) such that y(Tj) = Dj
C and y is a local minimum of the functional
C
C       E(y) = I[ K**2 ds],
C
C where I[_] is the integral with respect to arc length s,
C and K(s) denotes signed curvature.  The natural boundary
C conditions associated with minimizing E over nonlinear
C interpolating splines imply continuous curvature at the
C interior knots, K = 0 at the endpoints, and the Euler
C equation
C
C       K''(s) + 0.5*K**3 = 0
C
C is satisfied in each open interval (Ti,Ti+1).
C
C   The spline curve y models a thin beam (elastica) of
C homogeneous material and uniform cross sectional area
C with frictionlessly rotating sliders at the knots.  E(y)
C is proportional to the elastic strain energy.
C
C Reference:  Michael A. Malcolm, "On the computation of
C             nonlinear spline functions", SIAM J. Numerical
C             Analysis, Vol. 14, No. 2, April 1977.
C
C
C DISCRETIZATION:
C
C   Using K(t) = y''(t)/[1+y'(t)**2]**(3/2) and s'(t) =
C [1+y'(t)**2]**(1/2), we obtain
C
C       E(y) = I[ y''(t)**2/[1+y'(t)**2]**(5/2) dt],
C
C where I[_] denotes the integral from T1 to Tn.  The domain
C [T1,Tn] is uniformly partitioned with constant mesh width
C H, resulting in m grid points T1+(i-1)*H, i = 1,...,m with
C corresponding function values y(i) = y(T1+(i-1)*H).  The
C data points are perturbed if necessary so that, for j = 1,
C ...,n, Tj = T1+(i-1)*H for some integer i = IND(j), and
C y(i) = y(IND(j)) = y(Tj) = Dj.
C
C   We use central difference approximations to derivatives,
C and approximate the integral by the composite trapezoidal
C rule to obtain the discretized functional
C
C       phi(y) = H*Sum[ D2y(i)*a(i) ],
C
C where
C
C       D1y(i) = (y(i+1)-y(i-1))/(2*H),
C
C       D2y(i) = (y(i+1)-2*y(i)+y(i-1))/H**2,
C
C       a(i) = D2y(i)/[1+D1y(i)**2]**(5/2),
C
C and Sum[_] denotes the sum over i = 2 to m-1.
C
C   Denote by S0 the set of m-vectors h such that h(i) = 0
C if i = IND(j) for j = 1 to n.  This is the space of
C perturbations for y, chosen to preserve the interpolation
C conditions provided the initial estimate of y satisfies
C the conditions.  Let P0 denote the orthogonal projection
C onto S0.  P0 is applied by simply zeroing the appropriate
C components.
C
C The components of the gradient are
C
C grad(phi)(i) = (2/H)*[a(i-1)-2*a(i)+a(i+1)-b(i-1)+b(i+1)],
C
C where
C
C       b(i) = (5/4)*H*D1y(i)*D2y(i)*a(i)/[1+D1y(i)**2],
C
C for i = 2 to m-1.  grad(phi)(1) = grad(phi)(m) = 0 since
C the endpoint values y(1) and y(m) are fixed by the inter-
C polation conditions.  The gradient is projected onto S0 by
C zeroing the appropriate components.
C
C
C METHOD:
C
C   The functional phi is minimized subject to the interpo-
C lation constraints by the Polak-Ribiere nonlinear conjugate
C gradient method using the Sobolev gradient defined below.
C
C   The Sobolev inner product associated with a curve (the
C current approximation to the solution) y is
C
C       <g,h>_y = I[ g''(s)*h''(s) ds ],
C
C where s(t) is the arc length associated with y, and
C
C       g'(s) = g'(t)/(1+y'(t)**2)**(1/2),
C
C       g''(s) = [ (1+y'(t)**2)*g''(t) -
C                  y'(t)*y''(t)*g'(t) ]/(1+y'(t)**2)**2.
C
C Note that <g,h>_y is positive (defines an inner product)
C on functions with two or more zeros (and square integrable
C second derivatives).  The discretized inner product
C on S0 is
C
C       <g,h>_y = Sum[ Dg(i)*Dh(i) ],
C
C where the sum is over i = 2 to m-1, and Dg is the dis-
C cretization of
C
C       g''(s)*Sqrt(s') = [s'(t)*g''(t)-s''(t)*g'(t)]/
C                         [s'(t)**(5/2)]
C
C              = [ [1+y'(t)**2]*g''(t)-y'(t)*y''(t)*g'(t) ]/
C                [ [1+y'(t)**2]**(7/4) ].
C
C Denote the discretized L2 inner product by
C
C       <r,s>_(m-2) = Sum[ r(i)*s(i) ].
C
C with the sum over i = 2 to m-1.
C
C   Then <g,h>_y = <Dg,Dh>_(m-2) = <D**T*Dg,h>_m =
C <D**T*Dg,P0*h>_m = <P0*D**T*Dg,h>_S0, where P0 is the
C orthogonal projection onto S0 (and the adjoint of D is
C P0*D**T).  Now let g be the Sobolev gradient (S-gradient)
C of phi at y.  Then by the Rietz Representation Theorem,
C
C       phi'(y)h = <grad(phi),h>_S0 = <g,h>_y
C                = <P0*D**T*Dg,h>_S0
C
C for all h in S0, so that the Sobolev gradient g is defined
C by
C
C       DTD*g = grad(phi),
C
C where DTD is the restriction of P0*D**T*D to S0.
C
C   Note that D**T*D is a pentadiagonal symmetric positive
C semi-definite matrix of order m, but we implicitly omit
C the rows and columns associated with constraints (zeros in
C grad(phi)) to obtain a positive definite matrix of order
C m-n.  More precisely, we omit the first and last rows and
C columns, and treat the other constraints by zeroing the
C row and column elements other than the diagonal (which is
C retained in order to preserve the scaling and condition
C number).  The resulting order-(m-2) system is solved by a
C direct method (using a U**T*U factorization for upper
C triangular matrix U).  DTD may be thought of as a precon-
C ditioner for the standard gradient.
C
C
C REQUIRED LIBRARY:
C
C   This program must be linked with a library providing the
C LAPACK subroutine DPBSV and the double precision BLAS.
      SUBROUTINE CSPLIN (LDA,N,T,D, ABD, S,IER)
      INTEGER LDA, N, IER
      DOUBLE PRECISION T(N), D(N), ABD(LDA,N), S(N)
C
C***********************************************************
C
C                                              From DNSPLIN1
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   08/17/01
C
C   Given a sequence of N data points (T(i),D(i)) with
C strictly increasing abscissae, this routine solves the
C symmetric positive definite tridiagonal linear system for
C the knot-derivative values S(i) defining the C2 Hermite
C cubic spline interpolant with natural end conditions.
C
C On input:
C
C       LDA = Row dimension of ABD.  LDA >= 2.
C
C       N = Number of data points.  N >= 2.
C
C       T,D = Arrays of length N containing the data points
C             (abscissae and ordinates, respectively, of the
C             knots) with strictly increasing abscissae.  T
C             is not tested for validity.
C
C The above parameters are not altered by this routine.
C
C       ABD = Array dimensioned LDA by M for M >= N.
C
C       S = Array of length at least N.
C
C On output:
C
C       ABD = Garbage.
C
C       S = Knot-derivative values of the cubic spline if
C           IER = 0.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if LDA or N is outside its valid
C                     range.
C
C Modules required by CSPLIN:  None
C
C***********************************************************
C
      DOUBLE PRECISION A1, A2, DD1, DD2
      INTEGER J
C
C Test for valid input parameters.
C
      IF (LDA .LT. 2  .OR.  N .LT. 2) GO TO 11
C
C Store the symmetric tridiagonal matrix in ABD and right
C   hand side in S.  The first row of ABD contains the off-
C   diagonal elements 1/h(i) = 1/(T(i+1)-T(i)) in column
C   i+1 for i = 1 to N-1.  The second row contains the
C   diagonal elements 2/h(i-1) + 2/h(i) in column i for
C   i = 1 to N, where 2/h(0) = 2/h(N) = 0.  The right hand
C   side is S(i) = 3*[dd(i-1)/h(i-1) + dd(i)/h(i)], where
C   dd(i) = (D(i+1)-D(i))/h(i) and dd(0) = dd(N) = 0.
C
      A2 = 1.D0/(T(2)-T(1))
      DD2 = (D(2)-D(1))*A2*A2
      ABD(2,1) = 2.D0*A2
      S(1) = 3.D0*DD2
      DO 1 J = 2,N-1
        A1 = A2
        A2 = 1.D0/(T(J+1)-T(J))
        DD1 = DD2
        DD2 = (D(J+1)-D(J))*A2*A2
        ABD(1,J) = A1
        ABD(2,J) = 2.D0*(A1+A2)
        S(J) = 3.D0*(DD1+DD2)
    1   CONTINUE
      ABD(1,N) = A2
      ABD(2,N) = 2.D0*A2
      S(N) = 3.D0*DD2
C
C Forward elimination:
C
      DO 2 J = 2,N
        A1 = ABD(1,J)/ABD(2,J-1)
        ABD(2,J) = ABD(2,J) - A1*ABD(1,J)
        S(J) = S(J) - A1*S(J-1)
    2   CONTINUE
C
C Back solve:
C
      S(N) = S(N)/ABD(2,N)
      DO 3 J = N-1,1,-1
        S(J) = (S(J) - ABD(1,J+1)*S(J+1))/ABD(2,J)
    3   CONTINUE
C
C No error encountered.
C
      IER = 0
      RETURN
C
C Invalid input parameter.
C
   11 IER = 1
      RETURN
      END
      DOUBLE PRECISION FUNCTION DSTORE (X)
      DOUBLE PRECISION X
C
C***********************************************************
C
C                                                 From GLRH2
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   02/25/96
C
C   This function forces its argument X to be stored in a
C memory location, thus providing a means of determining
C floating point number characteristics (such as the machine
C precision) when it is necessary to avoid computation in
C high precision registers.
C
C On input:
C
C       X = Double precision value to be stored.
C
C X is not altered by this function.
C
C On output:
C
C       DSTORE = Value of X after it has been stored and
C                possibly truncated or rounded to the double
C                precision word length.
C
C Modules required by DSTORE:  None
C
C***********************************************************
C
      DOUBLE PRECISION Y
      COMMON/STCOM/Y
      Y = X
      DSTORE = Y
      RETURN
      END
      DOUBLE PRECISION FUNCTION FMIN (AX,BX,TOL,F, IFLG )
      DOUBLE PRECISION AX, BX, TOL, F
      INTEGER IFLG
C
C***********************************************************
C
C                                                  From GLRH
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   11/23/95
C
C   This function returns an approximation to the point at
C which a real-valued function F attains a minimum on the
C interval (AX,BX).  It is based on the function by the same
C name by R. Brent (01/01/73), modified to use reverse
C communication.
C
C   The method used is a combination of golden section
C search and successive parabolic interpolation.  Conver-
C gence is never much slower than that for a Fibonacci
C search.  If F has a continuous second derivative which is
C positive at the minimum (which is not at AX or BX), then
C convergence is superlinear, and usually of the order of
C about 1.324....
C
C   The function F is never evaluated at two points closer
C together than EPS*Abs(FMIN) + (TOL/3), where EPS is
C approximately the square root of the relative machine
C precision.  If F is a unimodal function and the computed
C values of F are always unimodal when separated by at least
C EPS*Abs(X*) + (TOL/3), then FMIN approximates the abscissa
C of the global minimum of F on the interval (AX,BX) with
C an error less than 3*EPS*Abs(FMIN) + TOL.  If F is not
C unimodal, then FMIN may approximate a local, but perhaps
C non-global, minimum to the same accuracy.
C
C On input:
C
C       AX,BX = Endpoints of the initial interval -- the
C               interval over which F is to be minimized.
C
C       TOL = Desired length of the interval of uncertainty
C             of the final result.  TOL .GE. 0.
C
C       F = Function value F(FMIN) if IFLG .NE. 0, or unused
C           dummy parameter on the first call.
C
C The above parameters are not altered by this function.
C
C       IFLG = Reverse communication flag:
C              IFLG = 0 if this is the first call to FMIN
C                       for a given minimization problem.
C              IFLG > 0 if F contains a function value
C                       at the point FMIN returned by the
C                       previously call.  The value of
C                       IFLG must be the value returned
C                       by the previous call.
C
C On Output:
C
C       IFLG = Reverse communication flag:
C              IFLG = 0 if FMIN is the solution.
C              IFLG > 0 if FMIN contains a point at which
C                       F is to be evaluated.  FMIN must
C                       be called again.
C              IFLG < 0 if FMIN was called with an invalid
C                       value of IFLG (neither 0 nor the
C                       returned value).
C
C       FMIN = Approximation to the point at which F is
C              minimized (IFLG = 0), or point at which F
C              is to be evaluated (IFLG > 0).
C
C Module required by FMIN:  DSTORE
C
C Intrinsic functions called by FMIN:  ABS, SIGN, SQRT
C
C Reference:  Richard Brent, Algorithms for Minimization
C             Without Derivatives, Prentice-Hall, Inc.
C             (1973).
C
C***********************************************************
C
      DOUBLE PRECISION A, B, C, D, E, FU, FV, FW, FX, EPS,
     .                 P, Q, R, TOL1, TOL2, U, V, W, X, XM
      DOUBLE PRECISION DSTORE
      SAVE
C
C Test IFLG.
C
      IF (IFLG .EQ. 1) GO TO 2
      IF (IFLG .EQ. 2) GO TO 6
      IF (IFLG .NE. 0) THEN
        IFLG = -ABS(IFLG)
        RETURN
      ENDIF
C
C C is the squared inverse of the golden ratio.
C
      C = 0.5D0*(3.D0 - SQRT(5.D0))
C
C EPS is approximately the square root of the relative
C   machine precision.
C
      EPS = 1.D0
    1 EPS = EPS/2.D0
        TOL1 = DSTORE(1.D0 + EPS)
        IF (TOL1 .GT. 1.D0) GO TO 1
      EPS = SQRT(EPS)
C
C Initialization:
C
      A = AX
      B = BX
      V = A + C*(B - A)
      W = V
      X = V
      E = 0.D0
C
C Get F(X):  IFLG = 1.
C
      IFLG = 1
      FMIN = X
      RETURN
    2 FX = F
      FV = FX
      FW = FX
C
C Main loop:
C
    3 XM = 0.5D0*(A + B)
        TOL1 = EPS*ABS(X) + TOL/3.D0
        TOL2 = 2.D0*TOL1
C
C Test for termination.
C
        IF (ABS(X - XM) .LE. (TOL2 - 0.5D0*(B - A)))
     .    GO TO 7
C
C Test for golden-section necessary.
C
        IF (ABS(E) .LE. TOL1) GO TO 4
C
C Fit a parabola.
C
        R = (X - W)*(FX - FV)
        Q = (X - V)*(FX - FW)
        P = (X - V)*Q - (X - W)*R
        Q = 2.D0*(Q - R)
        IF (Q .GT. 0.D0) P = -P
        Q =  ABS(Q)
        R = E
        E = D
C
C Test for parabola acceptable.
C
        IF ( ABS(P) .GE. ABS(0.5D0*Q*R)  .OR.
     .       P .LE. Q*(A - X) .OR.  P .GE. Q*(B - X) )
     .    GO TO 4
C
C Take a parabolic interpolation step.
C
        D = P/Q
        U = X + D
C
C F must not be evaluated too close to AX or BX.
C
        IF ((U - A) .LT. TOL2) D = SIGN(TOL1, XM - X)
        IF ((B - U) .LT. TOL2) D = SIGN(TOL1, XM - X)
        GO TO 5
C
C Take a golden-section step.
C
    4   IF (X .GE. XM) E = A - X
        IF (X .LT. XM) E = B - X
        D = C*E
C
C F must not be evaluated too close to X.
C
    5   IF (ABS(D) .GE. TOL1) THEN
          U = X + D
        ELSE
          U = X + SIGN(TOL1, D)
        ENDIF
C
C Get F(U).  IFLG = 2.
C
        IFLG = 2
        FMIN = U
        RETURN
    6   FU = F
C
C Update A, B, V, W, and X.
C
        IF (FU .LE. FX) THEN
          IF (U .GE. X) A = X
          IF (U .LT. X) B = X
          V = W
          FV = FW
          W = X
          FW = FX
          X = U
          FX = FU
          GO TO 3
        ENDIF
C
        IF (U .LT. X) THEN
          A = U
        ELSE
          B = U
        ENDIF
        IF (FU .LE. FW  .OR.  W .EQ. X) THEN
          V = W
          FV = FW
          W = U
          FW = FU
          GO TO 3
        ENDIF
        IF (FU .LE. FV  .OR.  V .EQ. X  .OR.  V .EQ. W) THEN
          V = U
          FV = FU
        ENDIF
        GO TO 3
C
C Return the solution.
C
    7 IFLG = 0
      FMIN = X
      RETURN
      END
      SUBROUTINE GRADL2 (N,H,IND,M,D1Y,D2Y,A, WK, G,GNRM1,
     .                   GNRMM,IER)
      INTEGER N, IND(N), M, IER
      DOUBLE PRECISION H, D1Y(M), D2Y(M), A(M), WK(M), G(M),
     .                 GNRM1, GNRMM
C
C***********************************************************
C
C                                              From DNSPLIN1
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   06/18/02
C
C   This subroutine computes the L2 gradient grad(phi) of
C the energy functional phi associated with a nonlinear
C spline interpolant (Function PHI).
C
C On input:
C
C       N = Number of data points.  2 <= N <= M-1.
C
C       H = Mesh width.  H > 0.
C
C       IND = Array of length N containing a strictly
C             increasing sequence of indexes (1 to M) of
C             the data values to be interpolated, with
C             IND(1) = 1 and IND(N) = M.
C
C       M = Number of grid points.  M = IND(N) >= 6.
C
C       D1Y,D2Y,A = Arrays of length M containing difference
C                   approximations to derivatives and terms
C                   defining phi(y) computed by Function
C                   PHI.
C
C The above parameters are not altered by this routine.
C
C       WK = Work array of length at least M.
C
C       G = Array of length at least M.
C
C On output:
C
C       WK = Garbage.
C
C       G = Components of the L2-gradient grad(phi) (partial
C           derivatives of phi with respect to the function
C           values Y(i)) unless IER > 0.  G is an element of
C           S0:  G(i) = 0 for i = IND(j), j = 1 to N.
C
C       GNRM1 = Vector 1-norm of grad(phi):  sum of absolute
C               values of the components.
C
C       GNRMM = Vector max-norm of grad(phi):  maximum of
C               the absolute values of the components.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if N, H, IND, or M is outside its
C                      valid range on input.  Output para-
C                      meters are not altered in this case.
C
C Modules required by GRADL2:  None
C
C Intrinsic functions called by GRADL2:  ABS, MAX
C
C***********************************************************
C
      DOUBLE PRECISION C1, C2, GMAX, SUM
      INTEGER I, J
C
C Test for invalid input parameters.
C
      IF (N .LT. 2  .OR.  H .LE. 0.D0  .OR.  M .LT. 6  .OR.
     .    M-N .LT. 1  .OR.  IND(1) .NE. 1  .OR.
     .    IND(N) .NE. M) GO TO 11
      DO 1 I = 2,N
        IF (IND(I) .LE. IND(I-1)) GO TO 11
    1   CONTINUE
C
C Store constants.
C
      C1 = 1.25D0*H
      C2 = 2.D0/H
C
C L2 gradient:  G(i) = (2/H)*[a(i-1)-2*a(i)+a(i+1)
C                             -b(i-1)+b(i+1)]
C
C   where a(i) = D2y(i)/[1+D1y(i)**2]**(5/2),
C
C         b(i) = (5/4)*H*D1y(i)*D2y(i)**2/
C                     [1+D1y(i)**2]**(7/2),
C
C for i = 2 to m-1.
C
C Store b(i) in WK.
C
      DO 2 I = 1,M
        WK(I) = C1*D1Y(I)*D2Y(I)*A(I)/(1.D0+D1Y(I)**2)
    2   CONTINUE
C
C Store grad(phi) in G.
C
      DO 3 I = 2,M-1
        G(I) = C2*(A(I-1)-2.D0*A(I)+A(I+1)-WK(I-1)+WK(I+1))
    3   CONTINUE
C
C Project G onto S0 by zeroing the components associated
C   with data points.
C
      DO 4 J = 1,N
        I = IND(J)
        G(I) = 0.D0
    4   CONTINUE
C
C Accumulate GNRM1 and GNRMM in SUM and GMAX, respectively.
C
      SUM = 0.D0
      GMAX = 0.D0
      DO 5 I = 2,M-1
        SUM = SUM + ABS(G(I))
        GMAX = MAX( GMAX,ABS(G(I)) )
    5   CONTINUE
      GNRM1 = SUM
      GNRMM = GMAX
C
C No errors encountered.
C
      IER = 0
      RETURN
C
C Invalid input parameter.
C
   11 IER = 1
      RETURN
      END
      SUBROUTINE GRADS4 (N,H,IND,M,D1Y,D2Y,LDA, ABD,G,
     .                   WK, IER)
      INTEGER N, IND(N), M, LDA, IER
      DOUBLE PRECISION H, D1Y(M), D2Y(M), ABD(LDA,M), G(M),
     .                 WK(M,3)
C
C***********************************************************
C
C                                              From DNSPLIN1
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/12/02
C
C   This subroutine applies a 4-th order smoother DTD**(-1)
C to the L2 gradient grad(phi), returning a Sobolev gradi-
C ent G of the functional phi associated with a nonlinear
C spline interpolant (Function PHI).  The linear system is
C solved by LAPACK Subroutine DPBSV.
C
C   DTD is the restriction to S0 of P0*D**T*D, where D is
C the discrete differential operator corresponding to Dg =
C g''(s)*Sqrt(s').  The Sobolev gradient is therefore the
C gradient associated with the following discretized inner
C product on S0:
C
C       <g,h>_y = <Dg,Dh>_(m-2),
C
C               = Sum[ Dg(i)*Dh(i) ],
C
C where the sum is over i = 2 to m-1.  This is the discrete
C analog of
C
C       <g,h>_y = I[ g''(s)*h''(s) ds ],
C
C where s(t) is the arc length associated with y, and
C
C       g'(s) = g'(t)/(1+y'(t)**2)**(1/2),
C
C       g''(s) = [ (1+y'(t)**2)*g''(t) -
C                  y'(t)*y''(t)*g'(t) ]/(1+y'(t)**2)**2.
C
C On input:
C
C       N = Number of data points.  2 <= N <= M-1.
C
C       H = Mesh width.  H > 0.
C
C       IND = Array of length N containing a strictly
C             increasing sequence of indexes (1 to M) of
C             the data values to be interpolated, with
C             IND(1) = 1 and IND(N) = M.
C
C       M = Number of grid points.  M = IND(N) >= 6.
C
C       D1Y,D2Y = Arrays of length M containing difference
C                 approximations to derivatives defining
C                 phi(y) (computed by Function PHI).
C
C       LDA = Row dimension of ABD.  LDA >= 3.
C
C The above parameters are not altered by this routine.
C
C       ABD = Array dimensioned LDA by M used to store the
C             symmetric positive semi-definite pentadiagonal
C             order-(M-2) matrix D**T*D.
C
C       G = Array of length M containing the components of
C           the L2-gradient grad(phi) (partial derivatives
C           of phi with respect to the function values
C           Y(i)).  G must be an element of S0:  G(i) = 0
C           for i = IND(j), j = 1 to N.  G may be computed
C           by Subroutine GRADL2.
C
C       WK = Work array of length at least 3*M.
C
C On output:
C
C       ABD = U**T*U factorization of the upper triangle
C             of D**T*D in LAPACK band storage format
C             (incomplete if IER = 2).
C
C       G = Sobolev gradient (an element of S0) unless
C           IER > 0, in which case G is not altered.
C
C       WK = Garbage.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if N, H, IND, M, or LDA is outside
C                     its valid range on input.  Output
C                     parameters are not altered in this
C                     case.
C             IER = 2 if DTD is not positive definite
C                     (singular to working precision).
C
C LAPACK subprograms required by GRADS4:  DPBSV, DPBTF2,
C                    DPBTRF, DPBTRS, DPOTF2, LSAME, XERBLA
C
C BLAS subprograms required by GRADS4:  DGEMM, DGEMV, DSYR,
C                                       DSYRK, DTBSV, DTRSM
C
C Intrinsic function called by GRADS4:  SQRT
C
C***********************************************************
C
      DOUBLE PRECISION H4, T1, T2
      INTEGER I, J
C
C Test for invalid input parameters.
C
      IF (N .LT. 2  .OR.  H .LE. 0.D0  .OR.  M .LT. 6  .OR.
     .    M-N .LT. 1  .OR.  IND(1) .NE. 1  .OR.
     .    IND(N) .NE. M  .OR.  LDA .LT. 3) GO TO 11
      DO 1 I = 2,N
        IF (IND(I) .LE. IND(I-1)) GO TO 11
    1   CONTINUE
C
C Compute terms defining D**T*D, where D is the second
C   derivative with respect to arc length scaled by
C   Sqrt(s')):
C
C   H4 = H**4
C
C   T1 = 1 + D1y(j)**2
C   T2 = 1/(H**4*(1+D1y(j)**2)**(7/2)) = 1/(H4*T1**3.5)
C
C   WK(j,1) = T1**2*T2
C   WK(j,2) = -H*T1*D1y(j)*D2y(j)*T2
C   WK(j,3) = (H*D1y(j)*D2y(j)/2)**2*T2
C
      H4 = H**4
      DO 2 J = 2,M-1
        T1 = 1.D0 + D1Y(J)*D1Y(J)
        T2 = 1.D0/(H4*T1**3.5D0)
        WK(J,1) = T1*T1*T2
        WK(J,2) = -H*T1*D1Y(J)*D2Y(J)*T2
        WK(J,3) = (H*D1Y(J)*D2Y(J)/2.D0)**2*T2
    2   CONTINUE
C
C Store the upper triangle of D**T*D in the first three rows
C   of ABD in LAPACK band storage format:
C
C       ABD(1,j+1) = (D**T*D)(j,j+2)  for j = 2 to m-3,
C       ABD(2,j) =   (D**T*D)(j,j+1)  for j = 2 to m-2,
C       ABD(3,j-1) = (D**T*D)(j,j)    for j = 2 to m-1.
C
      DO 3 J = 2,M-3
        ABD(1,J+1) = WK(J+1,1)-WK(J+1,3)
    3   CONTINUE
C
      DO 4 J = 2,M-2
        ABD(2,J) = -2.D0*WK(J,1)-WK(J,2)-2.D0*WK(J+1,1)+
     .             WK(J+1,2)
    4   CONTINUE
C
      ABD(3,1) = 4.D0*WK(2,1)+WK(3,1)-WK(3,2)+WK(3,3)
      DO 5 J = 3,M-2
        ABD(3,J-1) = WK(J-1,1)+WK(J-1,2)+WK(J-1,3)+
     .               4.D0*WK(J,1)+WK(J+1,1)-WK(J+1,2)+
     .               WK(J+1,3)
    5   CONTINUE
      ABD(3,M-2) = WK(M-2,1)+WK(M-2,2)+WK(M-2,3)+
     .             4.D0*WK(M-1,1)
C
C Convert D**T*D to DTD as follows:  the rows and columns
C   associated with a data point (zero element in G) are
C   zeroed except for the diagonal element which is retained
C   in order to preserve the scaling (and condition number).
C
      DO 6 J = 2,N-1
        I = IND(J)
        ABD(1,I-1) = 0.D0
        ABD(2,I-1) = 0.D0
        ABD(2,I) = 0.D0
        ABD(1,I+1) = 0.D0
    6   CONTINUE
C
C Overwrite ABD with a U**T*U factorization (for upper
C   triangular matrix U) stored in band format, and solve
C   the system.
C
      CALL DPBSV ('U',M-2,2,1,ABD,LDA,G(2),M-2,IER)
      IF (IER .NE. 0) THEN
        IER = 2
      ENDIF
      RETURN
C
C Invalid input parameter.
C
   11 IER = 1
      RETURN
      END
      SUBROUTINE LNSRCH (M,H,Y0,PHI0,P,FMTOL,OPT, S,Y,D1Y,
     .                   D2Y,A, PHIY,NEVAL,YMAX,DY,IER)
      INTEGER M, NEVAL, IER
      LOGICAL OPT
      DOUBLE PRECISION H, Y0(M), PHI0, P(M), FMTOL, S, Y(M),
     .                 D1Y(M), D2Y(M), A(M), PHIY, YMAX, DY
C
C***********************************************************
C
C                                              From DNSPLIN1
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   06/18/02
C
C   This subroutine minimizes phi(Y0+S*P) over positive
C step-sizes S (unless OPT = FALSE), where phi is the
C functional defined by Function PHI.
C
C On input:
C
C       M = Number of grid points.  M >= 2.
C
C       H = Mesh width.  H > 0.
C
C       Y0 = Array of length M containing the current
C            estimate of the minimizer of phi:  gridpoint
C            function values y(i).
C
C       PHI0 = phi(Y0).
C
C       P = Array of length M defining the search
C           direction for the descent method.
C
C       FMTOL = Nonnegative tolerance for Function FMIN:
C               the desired length of the interval of
C               uncertainty for the optimal step-size.
C
C       OPT = Flag with value TRUE iff the optimal step-
C             size is to be computed.
C
C The above parameters are not altered by this routine.
C
C       S = Step-size to be used if OPT = FALSE, or initial
C           estimate of the optimal step-size if OPT = TRUE.
C           The value S = 1 or the output from a previous
C           call is a reasonable choice for the initial
C           estimate.  S > 0.
C
C       Y = Array of length M.
C
C       D1Y,D2Y,A = Arrays of length at least M.
C
C On output:
C
C       S = Estimate of the optimal step-size (not altered
C           if OPT = FALSE).
C
C       Y = Minimizer of phi:  Y0 + S*P.
C
C       D1Y,D2Y,A = Finite difference approximations to
C                   derivatives and terms defining the
C                   integrand associated with phi(Y).
C
C       PHIY = phi(Y).
C
C       NEVAL = Number of evaluations of phi.
C
C       YMAX = Max-norm of Y.
C
C       DY = Relative change in Y:  Max-norm(Y-Y0)/(1+YMAX)
C            = S*Max-norm(P)/(1+YMAX).
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if M, H, FMTOL, or S is outside its
C                     valid range on input.
C             IER = 2 if P = 0.
C             Output parameters are not altered if IER > 0.
C
C Modules required by LNSRCH:  FMIN, PHI
C
C Intrinsic function called by LNSRCH:  ABS
C
C***********************************************************
C
      DOUBLE PRECISION PH2, PHS, S1, S2
      DOUBLE PRECISION FMIN, PHI
      INTEGER I, IFLG, NEV
C
C Test for invalid input.
C
      IF (M .LT. 2  .OR.  H .LE. 0.D0  .OR.
     .    FMTOL .LT. 0.D0  .OR.  S .LE. 0.D0) THEN
        NEVAL = 0
        IER = 1
        RETURN
      ENDIF
C
C Compute DY = Max-norm(P), and test for P = 0.
C
      DY = 0.D0
      DO 1 I = 1,M
        IF (ABS(P(I)) .GT. DY) DY = ABS(P(I))
    1   CONTINUE
      IF (DY .EQ. 0.D0) THEN
        NEVAL = 0
        IER = 2
        RETURN
      ENDIF
      NEV = 0
      IF (.NOT. OPT) GO TO 6
C
C Find a bracketing interval [S1,S2] = [0,S2] such that
C   phi(Y0+S2*P) > PHI0 = phi(Y0).  S2 is initialized to
C   2*S and doubled at each step as necessary.
C
      S1 = 0.D0
      S2 = S
    2 S2 = 2.D0*S2
        DO 3 I = 1,M
          Y(I) = Y0(I) + S2*P(I)
    3     CONTINUE
        PH2 = PHI (M,H,Y, D1Y,D2Y,A)
        NEV = NEV + 1
        IF (PH2 .LE. PHI0) GO TO 2
C
C Compute the optimal step-size S.
C
      IFLG = 0
      PHS = 0.D0
    4 S = FMIN (S1,S2,FMTOL,PHS, IFLG )
        IF (IFLG .EQ. 0) GO TO 6
        DO 5 I = 1,M
          Y(I) = Y0(I) + S*P(I)
    5     CONTINUE
        PHS = PHI (M,H,Y, D1Y,D2Y,A)
        NEV = NEV + 1
        GO TO 4
C
C Update Y and compute the final value PHIY = phi(Y), YMAX,
C   and DY = Max-norm(Y-Y0)/(1+YMAX) for Y = Y0+S*P.
C
    6 YMAX = 0.D0
      DO 7 I = 1,M
        Y(I) = Y0(I) + S*P(I)
        IF (ABS(Y(I)) .GT. YMAX) YMAX = ABS(Y(I))
    7   CONTINUE
      PHIY = PHI (M,H,Y, D1Y,D2Y,A)
      DY = S*DY/(1.D0+YMAX)
      NEVAL = NEV + 1
      IER = 0
      RETURN
      END
      DOUBLE PRECISION FUNCTION PHI (M,H,Y, D1Y,D2Y,A)
      INTEGER M
      DOUBLE PRECISION H, Y(M), D1Y(M), D2Y(M), A(M)
C
C***********************************************************
C
C                                              From DNSPLIN1
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   06/18/02
C
C   Given a discrete sequence of function values Y = y(1),
C y(2), ..., y(M) on a uniform grid with mesh width H, this
C function returns the value at Y of the discretized energy
C functional which approximates the integral with respect to
C arc length of squared curvature, along with the terms
C defining the functional and difference approximations to
C derivatives at the grid points:
C
C       PHI = H*Sum[ D2y(i)*a(i) ],
C
C where
C
C       D1y(i) = (y(i+1)-y(i-1))/(2*H),
C
C       D2y(i) = (y(i+1)-2*y(i)+y(i-1))/H**2,
C
C       a(i) = D2y(i)/[1+D1y(i)**2]**(5/2),
C
C and Sum[_] denotes the sum over i = 2 to m-1.  The end
C conditions give D2y(1) = D2y(m) = 0.
C
C On input:
C
C       M = Number of grid points.  M >= 2.
C
C       H = Mesh width.  H > 0.
C
C       Y = Array of length M containing the function
C           values y(i).
C
C The above parameters are not altered by this function.
C
C       D1Y,D2Y,A = Arrays of length at least M.
C
C On output:
C
C       D1Y,D2Y,A = Arrays containing the finite difference
C                   approximations defined above, along with
C                   the endpoint values associated with
C                   natural end conditions:  D2y(1) = D2y(M)
C                   = 0, and thus D1y(1) = (y(2)-y(1))/H and
C                   D1y(M) = (y(M)-y(M-1))/H.
C
C       PHI = Value of the discretized energy functional.
C
C The program is terminated with an error message written to
C the standard output unit if an input parameter is invalid.
C
C Modules required by PHI:  None
C
C Intrinsic function called by PHI:  SQRT
C
C***********************************************************
C
      DOUBLE PRECISION S1, S2, SUM, T
      INTEGER I
      IF (M .LT. 2  .OR.  H .LE. 0.D0) GO TO 11
C
C Initialization:
C
      S1 = 1.D0/(2.D0*H)
      S2 = 1.D0/(H*H)
C
C i = 1:  the terms defining phi(Y) are accumulated in SUM.
C
      D1Y(1) = (Y(2)-Y(1))/H
      D2Y(1) = 0.D0
      A(1) = 0.D0
      SUM = 0.D0
C
C i = 2 to m-1:
C
      DO 1 I = 2,M-1
        D1Y(I) = (Y(I+1)-Y(I-1))*S1
        D2Y(I) = (Y(I+1)-2.D0*Y(I)+Y(I-1))*S2
        T = 1.D0 + D1Y(I)*D1Y(I)
        A(I) = D2Y(I)/(T*T*SQRT(T))
        SUM = SUM + D2Y(I)*A(I)
    1   CONTINUE
C
C i = m:
C
      D1Y(M) = (Y(M)-Y(M-1))/H
      D2Y(M) = 0.D0
      A(M) = 0.D0
C
      PHI = H*SUM
      RETURN
C
C M < 2 or H <= 0.
C
   11 WRITE (*,100) M, H
  100 FORMAT (//' *** Error in PHI:  M = ',I4,', H = ',
     .        D10.3)
      STOP
      END
      SUBROUTINE PLTCRV (LUN,N,A,B,XS,YS,M,Y,PLT0,Y0,
     .                   PLTSIZ, IER)
      INTEGER LUN, N, M, IER
      LOGICAL PLT0
      DOUBLE PRECISION A, B, XS(N), YS(N), Y(M), Y0(M),
     .                 PLTSIZ
C
C***********************************************************
C
C                                              From DNSPLIN1
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   06/20/02
C
C   Given one or two sequences of function values y(1), ...,
C y(M), and optionally y0(1), ..., y0(M), associated with
C uniformly distributed grid points T(i) = A + (i-1)*H for
C H = (B-A)/(M-1), i = 1 to M, and a sequence of N data
C points (xs(i),ys(i)), this subroutine creates a level-2
C Encapsulated PostScript (EPS) file containing a plot of
C the curve (or curves) with marker symbols centered at the
C data points.  If the second curve y0 is included, it it
C drawn with dashed lines.  It is assumed that the data
C points lie inside the containing rectangle of the curves.
C
C On input:
C
C       LUN = Logical unit number in the range 0 to 99.
C             The unit should be opened with an appropriate
C             file name before the call to this routine.
C
C       N = Number of data points to be marked.
C           0 <= N <= M.
C
C       A,B = Left and right endpoints defining the domain.
C             B-A > 0.
C
C       XS,YS = Arrays of length N containing the abscissae
C               and ordinates of the data points to be
C               marked.
C
C       M = Number of grid points and function values.
C           M >= 2.
C
C       Y = Array of length M containing the function
C           values.
C
C       PLT0 = Logical variable with value TRUE iff the
C              second curve is to be drawn.
C
C       Y0 = Array of length M containing function values
C            defining the second curve if PLT0 = TRUE, or
C            dummy variable otherwise.
C
C       PLTSIZ = Plot size in inches.  The view volume
C                (bounding box defined by the function
C                domain and range) is mapped to a rectangu-
C                lar viewport with maximum side-length equal
C                to .88*PLTSIZ.  The containing rectangle is
C                centered on the 8.5 by 11 inch page, and
C                its boundary is drawn.  Labels below, to
C                the left of, and to the right of the rec-
C                tangle extend its dimensions by 1/4 inch,
C                2/3 inches, and 1/6 inch, respectively, for
C                FSIZ = 12 pts.  1.0 <= PLTSIZ <= 6.5.
C
C Input parameters are not altered by this routine.
C
C On output:
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if LUN, N, A, B, M, or PLTSIZ is
C                     outside its valid range.
C             IER = 2 if all function values are identical.
C             IER = 3 if an error was encountered in writing
C                     to unit LUN.
C
C   Various plotting options can be controlled by altering
C the data statement below.
C
C Modules required by PLTCRV:  None
C
C Intrinsic functions called by PLTCRV:  CHAR, MAX, MIN,
C                                        NINT, REAL
C
C***********************************************************
C
      INTEGER I, IH, IPX1, IPX2, IPY1, IPY2, IW
      LOGICAL ASPECT, AXES, SQUARE
      REAL    BSIZ, DASHL, DX, DY, FSIZ, H, R, SFX, SFY,
     .        T, TX, TY, WX1, WX2, WY1, WY2, X, XM, YM
C
      DATA    ASPECT/.TRUE./, AXES/.TRUE./, BSIZ/6.0/,
     .        DASHL/4.0/, FSIZ/12.0/, SQUARE/.FALSE./
C
C Local parameters:
C
C ASPECT =    Logical variable with value TRUE if the
C               aspect ratio of the view volume is to be
C               preserved in the viewport, FALSE if the
C               viewport is to be square.  Note the limit-
C               ation on the size of the aspect ratio R
C AXES =      Logical variable with value TRUE if the x and
C               y axes (intersected with the view volume)
C               are to be drawn
C BSIZ =      Box size in points for marker symbols
C DASHL =     Length (in points, at 72 points per inch) of
C               dashes and spaces in a dashed line pattern
C               used for connecting pseudo-vertices
C DX =        View volume width WX2-WX1
C DY =        View volume height WY2-WY1
C FSIZ =      Font size in points for labels
C H =         Mesh width
C I =         Array index
C IH =        Height of the viewport in points
C IPX1,IPY1 = X and y coordinates (in points) of the lower
C               left corner of the containing rectangle
C IPX2,IPY2 = X and y coordinates (in points) of the upper
C               right corner of the containing rectangle
C IW =        Width of the viewport in points
C R =         Aspect ratio of the viewport:
C               Max(1/2,Min(2,DX/DY)) if ASPECT = TRUE,
C               1 otherwise
C SFX,SFY =   Scale factors for mapping world coordinates
C               (window coordinates in [WX1,WX2] X [WY1,WY2])
C               to viewport coordinates
C SQUARE =    Logical variable with value TRUE if markers
C               are to be square boxes, FALSE if markers are
C               to be X's
C T =         Temporary variable
C TX,TY =     Translation vector for mapping world coordi-
C               nates to viewport coordinates
C WX1,WY1 =   X and y (world) coordinates of the lower left
C               corner of the window (view volume)
C WX2,WY2 =   X and y (world) coordinates of the upper right
C               corner of the window (view volume)
C X =         Grid point abscissae value
C XM,YM =     Location in world coordinates of the center of
C               a marker symbol or beginning of a label
C
C
C Test for error 1.
C
      IF (LUN .LT. 0  .OR.  LUN .GT. 99  .OR.  N .LT. 0
     .    .OR.  N .GT. M  .OR.  B-A .LE. 0.0  .OR.  M .LT. 2
     .    .OR.  PLTSIZ .LT. 1.0  .OR.  PLTSIZ .GT. 6.5D0)
     .   GO TO 11
C
C Compute the window (view volume) corner coordinates
C   (WX1,WY1) and WX2,WY2).
C
      WX1 = A
      WX2 = B
      WY1 = REAL(Y(1))
      WY2 = WY1
      DO 1 I = 2,M
        WY1 = MIN(WY1,REAL(Y(I)))
        WY2 = MAX(WY2,REAL(Y(I)))
    1   CONTINUE
      IF (PLT0) THEN
        DO 2 I = 1,M
          WY1 = MIN(WY1,REAL(Y0(I)))
          WY2 = MAX(WY2,REAL(Y0(I)))
    2     CONTINUE
      ENDIF
C
C Compute the dimensions and aspect ratio of the view
C   volume.
C
      DX = WX2 - WX1
      DY = WY2 - WY1
      IF (DY .LE. 0.0) GO TO 12
      R = 1.0
      IF (ASPECT) R = MAX(0.5,MIN(2.0,DX/DY))
C
C Compute the lower left (IPX1,IPY1) and upper right
C   (IPX2,IPY2) corner coordinates of the containing
C   rectangle.  The coordinates, specified in default user
C   space units (points, at 72 points/inch with origin at
C   the lower left corner of the page), are chosen to have
C   aspect ratio R, and to center the plot on the 8.5 by 11
C   inch page.  The center of the page is (306,396), and
C   T = PLTSIZ/2 in points.
C
      T = 36.0*PLTSIZ
      IF (R .GE. 1.0) THEN
        IPX1 = 306 - NINT(T)
        IPX2 = 306 + NINT(T)
        IPY1 = 396 - NINT(T/R)
        IPY2 = 396 + NINT(T/R)
      ELSE
        IPX1 = 306 - NINT(T*R)
        IPX2 = 306 + NINT(T*R)
        IPY1 = 396 - NINT(T)
        IPY2 = 396 + NINT(T)
      ENDIF
C
C Output header comments.  The bounding box corner coordi-
C   nates are obtained by extending the containing rectangle
C   to include its boundary and the labels.
C
      WRITE (LUN,100,ERR=13) IPX1-NINT(4.0*FSIZ),
     .                       IPY1-NINT(1.5*FSIZ),
     .                       IPX2+NINT(FSIZ), IPY2+1
  100 FORMAT ('%!PS-Adobe-3.0 EPSF-3.0'/
     .        '%%BoundingBox:',4I4/
     .        '%%Title:  Nonlinear Spline'/
     .        '%%Creator:  DNSPLIN1'/
     .        '%%EndComments')
C
C Set the line thickness to 2 points, and draw the boundary
C   of the containing rectangle.
C
      T = 2.0
      WRITE (LUN,110,ERR=13) T
      WRITE (LUN,120,ERR=13) IPX1, IPY1
      WRITE (LUN,130,ERR=13) IPX1, IPY2
      WRITE (LUN,130,ERR=13) IPX2, IPY2
      WRITE (LUN,130,ERR=13) IPX2, IPY1
      WRITE (LUN,140,ERR=13)
      WRITE (LUN,150,ERR=13)
  110 FORMAT (F12.6,' setlinewidth')
  120 FORMAT (2I4,' moveto')
  130 FORMAT (2I4,' lineto')
  140 FORMAT ('closepath')
  150 FORMAT ('stroke')
C
C Set IW and IH to the width and height of a viewport
C   obtained by shrinking the containing rectangle by 12%
C   in each direction.
C
      IW = NINT(0.88*REAL(IPX2-IPX1))
      IH = NINT(0.88*REAL(IPY2-IPY1))
C
C Set up a viewport mapping.
C
      SFX = REAL(IW)/DX
      SFY = REAL(IH)/DY
      TX = REAL(306-IW/2) - SFX*WX1
      TY = REAL(396-IH/2) - SFY*WY1
C
C Set the line thickness to 1 point.
C
      T = 1.0
      WRITE (LUN,110,ERR=13) T
C
C Draw the first curve Y (by creating and painting a path).
C
      H = DX/REAL(M-1)
      X = A
      WRITE (LUN,160,ERR=13) SFX*X+TX, SFY*REAL(Y(1))+TY
      DO 3 I = 2,M
        X = X + H
        WRITE (LUN,170,ERR=13) SFX*X+TX, SFY*REAL(Y(I))+TY
    3   CONTINUE
      WRITE (LUN,150,ERR=13)
  160 FORMAT (2F12.6,' moveto')
  170 FORMAT (2F12.6,' lineto')
      IF (PLT0) THEN
C
C Set the dashed line attribute.
C
        WRITE (LUN,180,ERR=13) DASHL, DASHL
  180   FORMAT ('[',2F12.6,'] 0 setdash')
C
C Draw the second curve Y0.
C
        X = A
        WRITE (LUN,160,ERR=13) SFX*X+TX, SFY*REAL(Y0(1))+TY
        DO 4 I = 2,M
          X = X + H
          WRITE (LUN,170,ERR=13) SFX*X+TX, SFY*REAL(Y0(I))+TY
    4     CONTINUE
        WRITE (LUN,150,ERR=13)
      ENDIF
C
C Restore the solid line attribute.
C
      WRITE (LUN,190,ERR=13)
  190 FORMAT ('[] 0 setdash')
      IF (AXES) THEN
C
C Draw x and y axes (if contained in the view volume).
C
        IF (WX1 .LT. 0.0  .AND.  WX2 .GT. 0.0) THEN
          WRITE (LUN,160,ERR=13) TX, SFY*WY1+TY
          WRITE (LUN,170,ERR=13) TX, SFY*WY2+TY
        ENDIF
        IF (WY1 .LT. 0.0  .AND.  WY2 .GT. 0.0) THEN
          WRITE (LUN,160,ERR=13) SFX*WX1+TX, TY
          WRITE (LUN,170,ERR=13) SFX*WX2+TX, TY
        ENDIF
      ENDIF
C
C Draw marker symbols (squares or X's) at the data points.
C
C   T = BSIZ/2,
C   (XM,YM) = Location of symbol center.
C
      T = BSIZ/2.0
      DO 5 I = 1,N
        XM = SFX*REAL(XS(I))+TX
        YM = SFY*REAL(YS(I))+TY
        WRITE (LUN,160,ERR=13) XM-T, YM-T
        IF (SQUARE) THEN
          WRITE (LUN,170,ERR=13) XM+T, YM-T
          WRITE (LUN,170,ERR=13) XM+T, YM+T
          WRITE (LUN,170,ERR=13) XM-T, YM+T
          WRITE (LUN,140,ERR=13)
        ELSE
          WRITE (LUN,170,ERR=13) XM+T, YM+T
          WRITE (LUN,160,ERR=13) XM+T, YM-T
          WRITE (LUN,170,ERR=13) XM-T, YM+T
        ENDIF
    5   CONTINUE
C
C Paint the path.
C
      WRITE (LUN,150,ERR=13)
C
C Select a font and scale it.
C
      T = FSIZ
      WRITE (LUN,200,ERR=13) T
  200 FORMAT ('/Helvetica findfont'/
     .        F12.6,' scalefont setfont')
C
C Draw tick marks and labels at the extreme points of the
C   domain and range.  The tick mark length is T/2.
C
      XM = SFX*WX1+TX
      YM = REAL(IPY1)
      WRITE (LUN,160,ERR=13) XM, YM
      WRITE (LUN,170,ERR=13) XM, YM-T/2.0
      XM = XM - 1.5*T
      YM = YM - 1.5*T
      WRITE (LUN,160,ERR=13) XM, YM
      WRITE (LUN,210,ERR=13) WX1
  210 FORMAT ('(',F6.2,') show')
C
      XM = SFX*WX2+TX
      YM = REAL(IPY1)
      WRITE (LUN,160,ERR=13) XM, YM
      WRITE (LUN,170,ERR=13) XM, YM-T/2.0
      XM = XM - 1.5*T
      YM = YM - 1.5*T
      WRITE (LUN,160,ERR=13) XM, YM
      WRITE (LUN,210,ERR=13) WX2
C
      XM = REAL(IPX1)
      YM = SFY*WY1+TY
      WRITE (LUN,160,ERR=13) XM, YM
      WRITE (LUN,170,ERR=13) XM-T/2.0, YM
      XM = XM - 4.0*T
      YM = YM - 0.25*T
      WRITE (LUN,160,ERR=13) XM, YM
      WRITE (LUN,210,ERR=13) WY1
C
      XM = REAL(IPX1)
      YM = SFY*WY2+TY
      WRITE (LUN,160,ERR=13) XM, YM
      WRITE (LUN,170,ERR=13) XM-T/2.0, YM
      XM = XM - 4.0*T
      YM = YM - 0.25*T
      WRITE (LUN,160,ERR=13) XM, YM
      WRITE (LUN,210,ERR=13) WY2
C
C Paint the path and output the showpage command and
C   end-of-file indicator.
C
      WRITE (LUN,220,ERR=13)
  220 FORMAT ('stroke'/
     .        'showpage'/
     .        '%%EOF')
C
C HP's interpreters require a one-byte End-of-PostScript-Job
C   indicator (to eliminate a timeout error message):
C   ASCII 4.
C
      WRITE (LUN,230,ERR=13) CHAR(4)
  230 FORMAT (A1)
C
C No error encountered.
C
      IER = 0
      RETURN
C
C Invalid input parameter LUN, N, A, B, M, or PLTSIZ.
C
   11 IER = 1
      RETURN
C
C DY = 0:  constant function.
C
   12 IER = 2
      RETURN
C
C Error writing to unit LUN.
C
   13 IER = 3
      RETURN
      END
      SUBROUTINE READF (LIN,NMAX, N,T,D,H,IND,IER)
      INTEGER LIN, NMAX, N, IND(NMAX), IER
      DOUBLE PRECISION T(NMAX), D(NMAX), H
C
C***********************************************************
C
C                                              From DNSPLIN1
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   01/02/02
C
C   This routine reads a data set consisting of the number
C of data points N, the abscissae T, ordinates D, and mesh
C width H for the finite difference grid.  The indexes IND
C of the data points in the uniform grid are also computed.
C Format I4 is used for N, format D22.15 for the other
C values.
C
C On input:
C
C       LIN = Logical unit number for input.  0 .LE. LIN
C             .LE. 99.
C
C       NMAX = Maximum value of N.  NMAX >= 2.
C
C The above parameters are not altered by this routine.
C
C       T,D,IND = Arrays of length at least NMAX.
C
C On output:
C
C       N = Number of data points.  2 <= N <= NMAX.
C
C       T,D = Ordered sequence of data abscissae and ordi-
C             nates, respectively.  The abscissae must be
C             strictly increasing with separation between
C             adjacent values greater than H/2.
C
C       H = Positive mesh width for the uniform finite
C           difference grid.  The local discretization error
C           is O(H**2), and the data abscissae are perturbed
C           if necessary so that T(i)-T(1) is a multiple of
C           H for i = 1 to N.
C
C       IND = Indexes of the (perturbed) data points in the
C             grid:
C               IND(i) = NINT((T(i)-T(1))/H) + 1
C             for i = 1 to N.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if LIN or NMAX is outside its
C                     valid range.
C             IER = 2 if a read error occurred.
C             IER = 3 if N < 2 or N > NMAX.
C             IER = 4 if H <= 0.
C             IER = 5 if the abscissae are not strictly
C                     increasing or the separation between
C                     adjacent abscissae is too small to
C                     produce distinct indexes.
C
C             The error conditions are tested in the above-
C             specified order, and processing is terminated
C             immediately on encountering an error.
C
C Modules required by READF:  None
C
C Intrinsic functions called by READF:  DBLE, NINT
C
C***********************************************************
C
      INTEGER I
C
C Input formats:
C
  100 FORMAT (I4)
  110 FORMAT (D22.15)
C
C Test for invalid input parameters.
C
      IF (LIN .LT. 0  .OR.  LIN .GT. 99  .OR.  NMAX .LT. 2)
     .  GO TO 11
C
C Read the data set and test for errors.
C
      READ (LIN,100,ERR=12) N
      IF (N .LT. 2  .OR.  N .GT. NMAX) GO TO 13
      READ (LIN,110,ERR=12) (T(I), I = 1,N)
      READ (LIN,110,ERR=12) (D(I), I = 1,N)
      READ (LIN,110,ERR=12) H
      IF (H .LE. 0.D0) GO TO 14
C
C Compute IND.
C
      IND(1) = 1
      DO 1 I = 2,N
        IND(I) = NINT((T(I)-T(1))/H) + 1
        IF (IND(I) .LT. IND(I-1)) GO TO 15
        T(I) = T(1) + DBLE(IND(I)-1)*H
    1   CONTINUE
C
C No error encountered.
C
      IER = 0
      RETURN
C
C Input parameter LIN or NMAX is outside its valid range.
C
   11 IER = 1
      RETURN
C
C Read error encountered.
C
   12 IER = 2
      RETURN
C
C N < 2 or N > NMAX.
C
   13 IER = 3
      RETURN
C
C H is not positive.
C
   14 IER = 4
      RETURN
C
C Abscissae not strictly increasing or too close together.
C
   15 IER = 5
      RETURN
      END
      SUBROUTINE WRITF (LOUT,M,Y,D1Y,D2Y)
      INTEGER LOUT, M
      DOUBLE PRECISION Y(M), D1Y(M), D2Y(M)
C
C***********************************************************
C
C                                              From DNSPLIN1
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   06/20/02
C
C   This routine creates an output data set consisting of
C the number of grid points m (format I5) followed by the
C sequence of gridpoint function values (y(i), i = 1,m),
C the first derivative values (D1y(i), i = 1,m), and the
C second derivative values (D2y(i), i = 1,m) written with
C format D22.15.
C
C On input:
C
C       LOUT = Logical unit number for writing the solution.
C              0 .LE. LOUT .LE. 99.
C
C       M = Number of grid points.  M > 0.
C
C       Y = Array of length M containing the function
C           values.
C
C       D1Y,D2Y = Arrays of length M containing difference
C                 approximations to derivatives defining
C                 phi(y) (computed by Function PHI).
C
C Input parameters are not altered by this routine.
C
C An error message is written to the standard output unit
C if an input parameter is invalid or a write error is
C encountered writing to unit LOUT.
C
C Modules required by WRITF:  None
C
C***********************************************************
C
      INTEGER I
C
C Output formats:
C
  100 FORMAT (I5)
  110 FORMAT (D22.15)
C
C Test for invalid parameters.
C
      IF (LOUT .LT. 0  .OR.  LOUT .GT. 99  .OR.  M .LE. 0)
     .  GO TO 11
C
C Create the data set.
C
      WRITE (LOUT,100,ERR=12) M
      DO 1 I = 1,M
        WRITE (LOUT,110,ERR=12) Y(I)
    1   CONTINUE
      DO 2 I = 1,M
        WRITE (LOUT,110,ERR=12) D1Y(I)
    2   CONTINUE
      DO 3 I = 1,M
        WRITE (LOUT,110,ERR=12) D2Y(I)
    3   CONTINUE
      RETURN
C
C Invalid input parameter.
C
   11 WRITE (*,210)
  210 FORMAT (///10X,'*** Error in WRITF:  invalid input ',
     .               'parameter ***')
      RETURN
C
C Error writing to unit LOUT.
C
   12 WRITE (*,220)
  220 FORMAT (///10X,'*** Error in WRITF writing to unit ',
     .               'LOUT ***')
      RETURN
      END
