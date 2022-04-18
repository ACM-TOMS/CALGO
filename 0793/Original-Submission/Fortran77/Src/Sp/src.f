c
c
      SUBROUTINE ABMOD(N,NCAPM,MCAP,EPS,IROUT,A,B,XIR,XII,IS,ALPHA,BETA,
     +                 NCAP,KOUNT,IERR,BE,X,W,E,P0,P1,P2)
c
c This routine generates recursion coefficients for the modified
c measure
c
c                     dlambda(t)/p(t),
c
c where
c
c         p(t) = product[(1+zeta(mu)*t)**is(mu)],
c
c the product being extended from  mu=1  to  mu=mcap. (If  mcap=0,
c the routine simply returns the recursion coefficients of the measure
c dlambda.) This is done by a discretization method based on Gauss
c quadrature rules relative to the measure  dlambda. It is assumed
c that all  zeta's  are either real or occurring in conjugate complex
c pairs which involve consecutive values of  mu.
c
c   Input:  n  - - -  the number of recursion coefficients desired;
c                     type integer
c           ncapm - - an upper bound on the number of terms in the
c                     discrete inner product used in the discretization
c                     method; normally the choice ncapm=400 should be
c                     satisfactory; type integer
c           mcap  - - the number of poles (not counting multiplicities);
c                     type integer
c           eps - - - a (relative) error tolerance; type real
c           irout - - an integer selecting the routine for generating
c                     the recursion coefficients of the discretized
c                     inner product; if  irout=1, the ORTHPOL routine
c                     sti  is selected, otherwise the routine  lancz
c           a,b - - - real arrays of dimension  ncapm  containing the
c                     recursion coefficients of the measure  dlambda
c           xir,xii - real arrays of dimension  mcap  containing the
c                     real and imaginary parts, respectively, of
c                     zeta(mu); if  zeta(mu)  is real, put  xii(mu)=0.
c           is  - - - an integer array of dimension  mcap  containing
c                     the respective multiplicities of the poles
c
c   Output: alpha,beta - arrays of dimension  n  containing the
c                     recursion coefficients of the modified measure
c           ncap  - - the number of terms in the discrete inner product
c                     that yields convergence; type integer
c           kount - - the number of iterations required in the
c                     discretization method; type integer
c           ierr  - - an error flag, where
c                     ierr=0 on normal return
c                     ierr=1 indicates an error condition in the
c                            ORTHPOL routine  gauss
c                     ierr=2 indicates an error condition in the
c                            ORTHPOL routine  sti  or  lancz, whichever
c                            is used
c                     ierr=ncapm if the method does not converge with
c                            ncap less than or equal to ncapm
c
c All remaining arrays in the calling sequence are internal working
c arrays.
c
c
c The common declaration below is to be used only if the parameter
c theta  in Fermi-Dirac and Bose-Einstein integrals is incorporated
c into a weight function as described at the end of Section 5 of the
c companion paper
c
c      common/par/theta
C     .. Scalar Arguments ..
      REAL EPS
      INTEGER IERR,IROUT,KOUNT,MCAP,N,NCAP,NCAPM
C     ..
C     .. Array Arguments ..  The upper bounds XII, XIR and IS are all
c                            MCAP, but this doesn't work in Fortran 77
c                            if MCAP = 0.
      REAL A(NCAPM),ALPHA(N),B(NCAPM),BE(N),BETA(N),E(NCAPM),P0(NCAPM),
     +     P1(NCAPM),P2(NCAPM),W(NCAPM),X(NCAPM),XII(*),XIR(*)
      INTEGER IS(*)
C     ..
C     .. Local Scalars ..
      REAL EPSM,P
      INTEGER IE,IERRG,IMU,INCAP,K,MU
C     ..
C     .. External Functions ..
      REAL R1MACH
      EXTERNAL R1MACH
C     ..
C     .. External Subroutines ..
      EXTERNAL GAUSS,LANCZ,STI
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
      integer i1mach, nout
      nout = i1mach(2)
C     ..
      IERR = 0
      EPSM = R1MACH(3)
      INCAP = 1
      DO 10 K = 1,N
          ALPHA(K) = A(K)
          BETA(K) = B(K)
   10 CONTINUE
      IF (MCAP.GT.0) THEN
          KOUNT = -1
          DO 20 K = 1,N
              BETA(K) = 0.
   20     CONTINUE
          NCAP = (2*N-1)/2
   30     DO 40 K = 1,N
              BE(K) = BETA(K)
   40     CONTINUE
          KOUNT = KOUNT + 1
          IF (KOUNT.GT.1) INCAP = 2** (KOUNT/5)*N
          NCAP = NCAP + INCAP
          IF (NCAP.GT.NCAPM) THEN
              IERR = NCAPM
              RETURN

          END IF
c
c Discretization of the inner product
c
          CALL GAUSS(NCAP,A,B,EPSM,X,W,IERRG,E)
          IF (IERRG.NE.0) THEN
              IERR = 1
              RETURN

          END IF

          DO 60 K = 1,NCAP
              P = 1.
              IMU = 0
              DO 50 MU = 1,MCAP
                  IF (IMU.EQ.0) THEN
                      IF (XII(MU).EQ.0.) THEN
                          P = ((1.+XIR(MU)*X(K))**IS(MU))*P

                      ELSE
                          P = (((1.+XIR(MU)*X(K))**2+
     +                        (XII(MU)*X(K))**2)**IS(MU))*P
                          IMU = 1
                      END IF

                  ELSE
                      IMU = 0
                  END IF

   50         CONTINUE
              W(K) = W(K)/P
c
c The following statement incorporates the parameter  theta
c into a weight function as described at the end of Section 5
c of the companion paper. It replaces the statement immediately
c preceding it.
c
c          w(k)=sqrt(1.+.5*theta*x(k))*w(k)/p
c
   60     CONTINUE
c
c Computation of the desired recursion coefficients
c
          IF (IROUT.EQ.1) THEN
              CALL STI(N,NCAP,X,W,ALPHA,BETA,IE,P0,P1,P2)
              IF (IE.NE.0) THEN
                  WRITE (nout,FMT=9000) IE
                  IERR = 2
                  RETURN

              END IF

          ELSE
              CALL LANCZ(N,NCAP,X,W,ALPHA,BETA,IE,P0,P1)
              IF (IE.NE.0) THEN
                  IERR = 2
                  RETURN

              END IF

          END IF

          DO 70 K = 1,N
              IF (ABS(BETA(K)-BE(K)).GT.EPS*BETA(K)) GO TO 30
   70     CONTINUE
      END IF

      RETURN

 9000 FORMAT (1X,'ie in sti=',I3)
      END
c
c
      SUBROUTINE GAUSS(N,ALPHA,BETA,EPS,ZERO,WEIGHT,IERR,E)
c
c Given  n  and a measure  dlambda, this routine generates the n-point
c Gaussian quadrature formula
c
c     integral over supp(dlambda) of f(x)dlambda(x)
c
c        = sum from k=1 to k=n of w(k)f(x(k)) + R(n;f).
c
c The nodes are returned as  zero(k)=x(k) and the weights as
c weight(k)=w(k), k=1,2,...,n. The user has to supply the recursion
c coefficients  alpha(k), beta(k), k=0,1,2,...,n-1, for the measure
c dlambda. The routine computes the nodes as eigenvalues, and the
c weights in term of the first component of the respective normalized
c eigenvectors of the n-th order Jacobi matrix associated with  dlambda.
c It uses a translation and adaptation of the algol procedure  imtql2,
c Numer. Math. 12, 1968, 377-383, by Martin and Wilkinson, as modified
c by Dubrulle, Numer. Math. 15, 1970, 450. See also Handbook for
c Autom. Comput., vol. 2 - Linear Algebra, pp.241-248, and the eispack
c routine  imtql2.
c
c        Input:  n - - the number of points in the Gaussian quadrature	
c                      formula; type integer
c                alpha,beta - - arrays of dimension  n  to be filled
c                      with the values of  alpha(k-1), beta(k-1), k=1,2,
c                      ...,n
c                eps - the relative accuracy desired in the nodes
c                      and weights
c
c        Output: zero- array of dimension  n  containing the Gaussian
c                      nodes (in increasing order)  zero(k)=x(k), k=1,2,
c                      ...,n
c                weight - array of dimension  n  containing the
c                      Gaussian weights  weight(k)=w(k), k=1,2,...,n
c                ierr- an error flag equal to  0  on normal return,
c                      equal to  i  if the QR algorithm does not
c                      converge within 30 iterations on evaluating the
c                      i-th eigenvalue, equal to  -1  if  n  is not in
c                      range, and equal to  -2  if one of the beta's is
c                      negative.
c
c The array  e  is needed for working space.
c
C     .. Scalar Arguments ..
      REAL EPS
      INTEGER IERR,N
C     ..
C     .. Array Arguments ..
      REAL ALPHA(N),BETA(N),E(N),WEIGHT(N),ZERO(N)
C     ..
C     .. Local Scalars ..
      REAL B,C,F,G,P,R,S
      INTEGER I,II,J,K,L,M,MML
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,SIGN,SQRT
C     ..
      IF (N.LT.1) THEN
          IERR = -1
          RETURN

      END IF

      IERR = 0
      ZERO(1) = ALPHA(1)
      IF (BETA(1).LT.0.) THEN
          IERR = -2
          RETURN

      END IF

      WEIGHT(1) = BETA(1)
      IF (N.EQ.1) RETURN
      WEIGHT(1) = 1.
      E(N) = 0.
      DO 10 K = 2,N
          ZERO(K) = ALPHA(K)
          IF (BETA(K).LT.0.) THEN
              IERR = -2
              RETURN

          END IF

          E(K-1) = SQRT(BETA(K))
          WEIGHT(K) = 0.
   10 CONTINUE
      DO 80 L = 1,N
          J = 0
c
c Look for a small subdiagonal element.
c
   20     DO 30 M = L,N
              IF (M.EQ.N) GO TO 40
              IF (ABS(E(M)).LE.EPS* (ABS(ZERO(M))+ABS(ZERO(M+
     +            1)))) GO TO 40
   30     CONTINUE
   40     P = ZERO(L)
          IF (M.EQ.L) GO TO 80
          IF (J.EQ.30) GO TO 120
          J = J + 1
c
c Form shift.
c
          G = (ZERO(L+1)-P)/ (2.*E(L))
          R = SQRT(G*G+1.)
          G = ZERO(M) - P + E(L)/ (G+SIGN(R,G))
          S = 1.
          C = 1.
          P = 0.
          MML = M - L
c
c For i=m-1 step -1 until l do ...
c
          DO 70 II = 1,MML
              I = M - II
              F = S*E(I)
              B = C*E(I)
              IF (ABS(F).LT.ABS(G)) GO TO 50
              C = G/F
              R = SQRT(C*C+1.)
              E(I+1) = F*R
              S = 1./R
              C = C*S
              GO TO 60

   50         S = F/G
              R = SQRT(S*S+1.)
              E(I+1) = G*R
              C = 1./R
              S = S*C
   60         G = ZERO(I+1) - P
              R = (ZERO(I)-G)*S + 2.*C*B
              P = S*R
              ZERO(I+1) = G + P
              G = C*R - B
c
c Form first component of vector.
c
              F = WEIGHT(I+1)
              WEIGHT(I+1) = S*WEIGHT(I) + C*F
              WEIGHT(I) = C*WEIGHT(I) - S*F
   70     CONTINUE
          ZERO(L) = ZERO(L) - P
          E(L) = G
          E(M) = 0.
          GO TO 20

   80 CONTINUE
c
c Order eigenvalues and eigenvectors.
c
      DO 100 II = 2,N
          I = II - 1
          K = I
          P = ZERO(I)
          DO 90 J = II,N
              IF (ZERO(J).GE.P) GO TO 90
              K = J
              P = ZERO(J)
   90     CONTINUE
          IF (K.EQ.I) GO TO 100
          ZERO(K) = ZERO(I)
          ZERO(I) = P
          P = WEIGHT(I)
          WEIGHT(I) = WEIGHT(K)
          WEIGHT(K) = P
  100 CONTINUE
      DO 110 K = 1,N
          WEIGHT(K) = BETA(1)*WEIGHT(K)*WEIGHT(K)
  110 CONTINUE
      RETURN
c
c Set error - no convergence to an eigenvalue after 30 iterations.
c
  120 IERR = L
      RETURN

      END
c
c
      SUBROUTINE GCHRIS(N,IOPT,A,B,X,HP,HN,ALPHA,BETA,IERR)
c
c This routine implements the generalized Christoffel theorem in two
c special cases. Given the first  n  recursion coefficients for some
c measure  dlambda(t), it generates the first  n  recursion coefficients
c for the measure
c
c                dlambda(t)/(t-x)        if iopt=1
c                dlambda(t)/(t**2-x**2)  if iopt=2
c
c where  x  is a real number outside the support interval of  dlambda.
c The case  iopt=1  is identical with option 4 of the ORTHPOL routine
c chri.
c
c On entry,
c
c    n  - - is the number of recursion coefficients desired;
c           type integer
c    iopt - is an integer selecting the desired modification
c           of the measure  dlambda
c    a,b  - are real arrays of dimension  n  containing the
c           first  n  recursion coefficients for the (monic)
c           orthogonal polynomials relative to the measure
c           dlambda
c    x  - - is a real parameter defining the linear resp.
c           quadratic modification of the measure  dlambda
c    hp - - is the value of the integral of dlambda(t)/(x-t);
c           type real
c    hn - - is the value of the integral of dlmabda(t)/(-x-t);
c           type real
c
c On retrun,
c
c    alpha,beta - are real arrays of dimension  n  containing
c           the first  n  recursion coefficients of the modified
c           measure
c    ierr - is an error flag, type integer, where
c           ierr=0 on normal return
c           ierr-1 if  iopt is not in range
c
C     .. Scalar Arguments ..
      REAL HN,HP,X
      INTEGER IERR,IOPT,N
C     ..
C     .. Array Arguments ..
      REAL A(N),ALPHA(N),B(N),BETA(N)
C     ..
C     .. Local Scalars ..
      REAL E,EH,Q,QH
      INTEGER K
C     ..
      IERR = 0
      IF (IOPT.EQ.1) THEN
          ALPHA(1) = X - B(1)/HP
          BETA(1) = -HP
          Q = -B(1)/HP
          DO 10 K = 2,N
              E = A(K-1) - X - Q
              BETA(K) = Q*E
              Q = B(K)/E
              ALPHA(K) = Q + E + X
   10     CONTINUE
          RETURN

      ELSE IF (IOPT.EQ.2) THEN
          ALPHA(1) = X* (HP+HN)/ (HP-HN)
          BETA(1) = - (HP-HN)/ (2.*X)
          Q = -B(1)/HP
          QH = -HP/BETA(1)
          E = 0.
          DO 20 K = 2,N
              EH = Q + E + 2.*X - QH
              BETA(K) = QH*EH
              E = A(K-1) - X - Q
              QH = Q*E/EH
              ALPHA(K) = QH + EH - X
              Q = B(K)/E
   20     CONTINUE
          RETURN

      ELSE
          IERR = 1
          RETURN

      END IF

      END
c
c
      SUBROUTINE GQRAT(N,MCAP,ALPHA,BETA,XIR,XII,IS,ZG,WG,CONST,IERR,E)
c
c Given the first  n+1  recursion coefficients  alpha,beta  of the
c modified measure  dlambda(t)/p(t), this routine generates the n-point
c Gauss-type rational quadrature rule for the measure  dlambda,
c including its error constant
c
c     .. Scalar Arguments ..
      REAL CONST
      INTEGER IERR,MCAP,N
c     ..
c     .. Array Arguments ..  The upper bounds XII, XIR and IS are all
c                            MCAP, but this doesn't work in Fortran 77
c                            if MCAP = 0.
      REAL ALPHA(N+1),BETA(N+1),E(N),WG(N),XII(*),XIR(*),ZG(N)
      INTEGER IS(*)
c     ..
c     .. Local Scalars ..
      REAL EPSM,P
      INTEGER IERRG,IMU,K,MU
c     ..
c     .. External Functions ..
      REAL R1MACH
      EXTERNAL R1MACH
c     ..
c     .. External Subroutines ..
      EXTERNAL GAUSS
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC REAL
      integer i1mach, nout
      nout = i1mach(2)
c     ..
      IERR = 0
      EPSM = R1MACH(3)
c
c compute the n-point Gaussian quadrature rule for the modified measure
c
      CALL GAUSS(N,ALPHA,BETA,EPSM,ZG,WG,IERRG,E)
      IF (IERRG.NE.0) THEN
          WRITE (nout,FMT=9000) IERRG
          IERR = 1
          RETURN

      END IF
c
c compute the error constant and the weights of the desired rational
c Gauss-type quadrature rule
c
      CONST = BETA(1)
      DO 20 K = 1,N
          CONST = BETA(K+1)*CONST/REAL(2*K* (2*K-1))
          IF (MCAP.GT.0) THEN
              P = 1.
              IMU = 0
              DO 10 MU = 1,MCAP
                  IF (IMU.EQ.0) THEN
                      IF (XII(MU).EQ.0.) THEN
                          P = ((1.+XIR(MU)*ZG(K))**IS(MU))*P

                      ELSE
                          P = (((1.+XIR(MU)*ZG(K))**2+
     +                        (XII(MU)*ZG(K))**2)**IS(MU))*P
                          IMU = 1
                      END IF

                  ELSE
                      IMU = 0
                  END IF

   10         CONTINUE
              WG(K) = P*WG(K)
          END IF

   20 CONTINUE
      RETURN

 9000 FORMAT (1X,'ierrg in gauss=',I5)
      END
c
c
      SUBROUTINE KNUM(N,NU0,NUMAX,Z,EPS,A,B,RHO,NU,IERR,ROLD)
c
c This routine generates
c
c   rho(k)(z)=integral pi(k)(t)dlambda(t)/(z-t), k=0,1,2,...,n,
c
c where  pi(k)(t)  is the (monic) k-th degree orthogonal polynomial
c with respect to the measure  dlambda(t), and the integral is extended
c over the support of  dlambda. It is assumed that  z  is a complex
c number outside the smallest interval containing the support of
c dlambda. The quantities  rho(k)(z)  are computed as the first  n+1
c members of the minimal solution of the basic three-term recurrence
c relation
c
c      y(k+1)(z)=(z-a(k))y(k)(z)-b(k)y(k-1)(z), k=0,1,2,...,
c
c satisfied by the orthogonal polynomials  pi(k)(z).
c
c   Input:  n  - -  the largest integer  k  for which  rho(k)  is
c                   desired
c           nu0  -  an estimate of the starting backward recurrence
c                   index; if no better estimate is known, set
c                   nu0 = 3*n/2; for Jacobi, Laguerre and Hermite
c                   weight functions, estimates of  nu0  are generated
c                   respectively by the routines  nu0jac,nu0lag  and
c                   nu0her
c           numax - an integer larger than  n  cutting off backward
c                   recursion in case of nonconvergence; if  nu0
c                   exceeds  numax, then the routine aborts with the
c                   error flag  ierr  set equal to  nu0
c           z - - - the variable in  rho(k)(z); type complex
c           eps - - the relative accuracy to which the  rho(k)  are
c                   desired
c           a,b - - arrays of dimension  numax  to be supplied with the
c                   recurrence coefficients  a(k-1), b(k-1), k=1,2,...,
c                   numax.
c
c   Output: rho - - an array of dimension  n+1  containing the results
c                   rho(k)=rho(k-1)(z), k=1,2,...,n+1; type complex
c           nu  - - the starting backward recurrence index that yields
c                   convergence
c           ierr  - an error flag equal to zero on normal return, equal
c                   to  nu0  if  nu0 > numax, and equal to  numax in
c                   case of nonconvergence.
c
c The complex array  rold  of dimension  n+1  is used for working space.
c
c
c The arrays  rho,rold  are assumed to have dimension  n+1.
c
C     .. Scalar Arguments ..
      COMPLEX Z
      REAL EPS
      INTEGER IERR,N,NU,NU0,NUMAX
C     ..
C     .. Array Arguments ..
      COMPLEX RHO(*),ROLD(*)
      REAL A(NUMAX),B(NUMAX)
C     ..
C     .. Local Scalars ..
      COMPLEX R
      INTEGER J,J1,K,NP1
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,CMPLX
C     ..
      IERR = 0
      NP1 = N + 1
      IF (NU0.GT.NUMAX) THEN
          IERR = NU0
          RETURN

      END IF

      IF (NU0.LT.NP1) NU0 = NP1
      NU = NU0 - 5
      DO 10 K = 1,NP1
          RHO(K) = (0.,0.)
   10 CONTINUE
   20 NU = NU + 5
      IF (NU.GT.NUMAX) THEN
          IERR = NUMAX
          GO TO 60

      END IF

      DO 30 K = 1,NP1
          ROLD(K) = RHO(K)
   30 CONTINUE
      R = (0.,0.)
      DO 40 J = 1,NU
          J1 = NU - J + 1
          R = CMPLX(B(J1),0.)/ (Z-CMPLX(A(J1),0.)-R)
          IF (J1.LE.NP1) RHO(J1) = R
   40 CONTINUE
      DO 50 K = 1,NP1
          IF (ABS(RHO(K)-ROLD(K)).GT.EPS*ABS(RHO(K))) GO TO 20
   50 CONTINUE
   60 IF (N.EQ.0) RETURN
      DO 70 K = 2,NP1
          RHO(K) = RHO(K)*RHO(K-1)
   70 CONTINUE
      RETURN

      END
c
c
      SUBROUTINE LANCZ(N,NCAP,X,W,ALPHA,BETA,IERR,P0,P1)
c
c This routine carries out the same task as the routine  sti, but
c uses the more stable Lanczos method. The meaning of the input
c and output parameters is the same as in the routine  sti. (This
c routine is adapted from the routine RKPW in W.B. Gragg and
c W.J. Harrod,The numerically stable reconstruction of Jacobi
c matrices from spectral data'', Numer. Math. 44, 1984, 317-335.)
c
C     .. Scalar Arguments ..
      INTEGER IERR,N,NCAP
C     ..
C     .. Array Arguments ..
      REAL ALPHA(N),BETA(N),P0(NCAP),P1(NCAP),W(NCAP),X(NCAP)
C     ..
C     .. Local Scalars ..
      REAL GAM,PI,RHO,SIG,T,TK,TMP,TSIG,XLAM
      INTEGER I,K
C     ..
      IF (N.LE.0 .OR. N.GT.NCAP) THEN
          IERR = 1
          RETURN

      ELSE
          IERR = 0
      END IF

      DO 10 I = 1,NCAP
          P0(I) = X(I)
          P1(I) = 0.
   10 CONTINUE
      P1(1) = W(1)
      DO 30 I = 1,NCAP - 1
          PI = W(I+1)
          GAM = 1.
          SIG = 0.
          T = 0.
          XLAM = X(I+1)
          DO 20 K = 1,I + 1
              RHO = P1(K) + PI
              TMP = GAM*RHO
              TSIG = SIG
              IF (RHO.LE.0.) THEN
                  GAM = 1.
                  SIG = 0.

              ELSE
                  GAM = P1(K)/RHO
                  SIG = PI/RHO
              END IF

              TK = SIG* (P0(K)-XLAM) - GAM*T
              P0(K) = P0(K) - (TK-T)
              T = TK
              IF (SIG.LE.0.) THEN
                  PI = TSIG*P1(K)

              ELSE
                  PI = (T**2)/SIG
              END IF

              TSIG = SIG
              P1(K) = TMP
   20     CONTINUE
   30 CONTINUE
      DO 40 K = 1,N
          ALPHA(K) = P0(K)
          BETA(K) = P1(K)
   40 CONTINUE
      RETURN

      END
c
c
      SUBROUTINE RECUR(N,IPOLY,AL,BE,A,B,IERR)
c
c This subroutine generates the coefficients  a(k),b(k), k=0,1,...,n-1,
c in the recurrence relation
c
c       p(k+1)(x)=(x-a(k))*p(k)(x)-b(k)*p(k-1)(x),
c                            k=0,1,...,n-1,
c
c       p(-1)(x)=0,  p(0)(x)=1,
c
c for some classical (monic) orthogonal polynomials, and sets  b(0)
c equal to the total mass of the weight distribution. The results are
c stored in the arrays  a,b,  which hold, respectively, the coefficients
c a(k-1),b(k-1), k=1,2,...,n.
c
c       Input:  n - - the number of recursion coefficients desired
c               ipoly-integer identifying the polynomial as follows:
c                     1=Legendre polynomial on (-1,1)
c                     2=Legendre polynomial on (0,1)
c                     3=Chebyshev polynomial of the first kind
c                     4=Chebyshev polynomial of the second kind
c                     5=Jacobi polynomial with parameters  al=-.5,be=.5
c                     6=Jacobi polynomial with parameters  al,be
c                     7=generalized Laguerre polynomial with
c                       parameter  al
c                     8=Hermite polynomial
c               al,be-input parameters for Jacobi and generalized
c                     Laguerre polynomials
c
c       Output: a,b - arrays containing, respectively, the recursion
c                     coefficients  a(k-1),b(k-1), k=1,2,...,n.
c               ierr -an error flag, equal to  0  on normal return,
c                     equal to  1  if  al  or  be  are out of range
c                     when  ipoly=6  or  ipoly=7, equal to  2  if  b(0)
c                     overflows when  ipoly=6  or  ipoly=7, equal to  3
c                     if  n  is out of range, and equal to  4  if  ipoly
c                     is not an admissible integer. In the case  ierr=2,
c                     the coefficient  b(0)  is set equal to the largest
c                     machine-representable number.
c
c The subroutine calls for the function subroutines  r1mach,gamma  and
c alga. The routines  gamma  and  alga , which are included in this
c file, evaluate respectively the gamma function and its logarithm for
c positive arguments. They are used only in the cases  ipoly=6  and
c ipoly=7.
c
C     .. Scalar Arguments ..
      REAL AL,BE
      INTEGER IERR,IPOLY,N
C     ..
C     .. Array Arguments ..
      REAL A(N),B(N)
C     ..
C     .. Local Scalars ..
      REAL AL2,ALMACH,ALPBE,BE2,FKM1,T
      INTEGER K
C     ..
C     .. External Functions ..
      REAL ALGA,GAMMA,R1MACH
      EXTERNAL ALGA,GAMMA,R1MACH
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ATAN,EXP,LOG,REAL,SQRT
C     ..
      IF (N.LT.1) THEN
          IERR = 3
          RETURN

      END IF

      ALMACH = LOG(R1MACH(2))
      IERR = 0
      DO 10 K = 1,N
          A(K) = 0.
   10 CONTINUE
      IF (IPOLY.EQ.1) THEN
          B(1) = 2.
          IF (N.EQ.1) RETURN
          DO 20 K = 2,N
              FKM1 = REAL(K-1)
              B(K) = 1./ (4.-1./ (FKM1*FKM1))
   20     CONTINUE
          RETURN

      ELSE IF (IPOLY.EQ.2) THEN
          A(1) = .5
          B(1) = 1.
          IF (N.EQ.1) RETURN
          DO 30 K = 2,N
              A(K) = .5
              FKM1 = REAL(K-1)
              B(K) = .25/ (4.-1./ (FKM1*FKM1))
   30     CONTINUE
          RETURN

      ELSE IF (IPOLY.EQ.3) THEN
          B(1) = 4.*ATAN(1.)
          IF (N.EQ.1) RETURN
          B(2) = .5
          IF (N.EQ.2) RETURN
          DO 40 K = 3,N
              B(K) = .25
   40     CONTINUE
          RETURN

      ELSE IF (IPOLY.EQ.4) THEN
          B(1) = 2.*ATAN(1.)
          IF (N.EQ.1) RETURN
          DO 50 K = 2,N
              B(K) = .25
   50     CONTINUE
          RETURN

      ELSE IF (IPOLY.EQ.5) THEN
          B(1) = 4.*ATAN(1.)
          A(1) = .5
          IF (N.EQ.1) RETURN
          DO 60 K = 2,N
              B(K) = .25
   60     CONTINUE
          RETURN

      ELSE IF (IPOLY.EQ.6) THEN
          IF (AL.LE.-1. .OR. BE.LE.-1.) THEN
              IERR = 1
              RETURN

          ELSE
              ALPBE = AL + BE
              A(1) = (BE-AL)/ (ALPBE+2.)
              T = (ALPBE+1.)*LOG(2.) + ALGA(AL+1.) + ALGA(BE+1.) -
     +            ALGA(ALPBE+2.)
              IF (T.GT.ALMACH) THEN
                  IERR = 2
                  B(1) = R1MACH(2)

              ELSE
                  B(1) = EXP(T)
              END IF

              IF (N.EQ.1) RETURN
              AL2 = AL*AL
              BE2 = BE*BE
              A(2) = (BE2-AL2)/ ((ALPBE+2.)* (ALPBE+4.))
              B(2) = 4.* (AL+1.)* (BE+1.)/ ((ALPBE+3.)* (ALPBE+2.)**2)
              IF (N.EQ.2) RETURN
              DO 70 K = 3,N
                  FKM1 = REAL(K-1)
                  A(K) = .25* (BE2-AL2)/ (FKM1*FKM1* (1.+.5*ALPBE/FKM1)*
     +                   (1.+.5* (ALPBE+2.)/FKM1))
                  B(K) = .25* (1.+AL/FKM1)* (1.+BE/FKM1)*
     +                   (1.+ALPBE/FKM1)/ ((1.+.5* (ALPBE+1.)/FKM1)*
     +                   (1.+.5* (ALPBE-1.)/FKM1)*
     +                   (1.+.5*ALPBE/FKM1)**2)
   70         CONTINUE
              RETURN

          END IF

      ELSE IF (IPOLY.EQ.7) THEN
          IF (AL.LE.-1.) THEN
              IERR = 1
              RETURN

          ELSE
              A(1) = AL + 1.
              B(1) = GAMMA(AL+1.,IERR)
              IF (IERR.EQ.2) B(1) = R1MACH(2)
              IF (N.EQ.1) RETURN
              DO 80 K = 2,N
                  FKM1 = REAL(K-1)
                  A(K) = 2.*FKM1 + AL + 1.
                  B(K) = FKM1* (FKM1+AL)
   80         CONTINUE
              RETURN

          END IF

      ELSE IF (IPOLY.EQ.8) THEN
          B(1) = SQRT(4.*ATAN(1.))
          IF (N.EQ.1) RETURN
          DO 90 K = 2,N
              B(K) = .5*REAL(K-1)
   90     CONTINUE
          RETURN

      ELSE
          IERR = 4
      END IF

      END

      REAL FUNCTION ALGA(X)
c
c This is an auxiliary function subroutine (not optimized in any
c sense) evaluating the logarithm of the gamma function for positive
c arguments  x. It is called by the subroutine  gamma. The integer  m0
c in the first executable statement is the smallest integer  m  such
c that  1*3*5* ... *(2*m+1)/(2**m)  is greater than or equal to the
c largest machine-representable number. The routine is based on a
c rational approximation valid on [.5,1.5] due to W.J. Cody and
c K.E. Hillstrom; see Math. Comp. 21, 1967, 198-203, in particular the
c case  n=7  in Table II. For the computation of  m0  it calls upon the
c function subroutines  t  and  r1mach. The former, appended below,
c evaluates the inverse function  t = t(y)  of  y = t ln t.
c
C     .. Scalar Arguments ..
      REAL X
C     ..
C     .. Local Scalars ..
      REAL P,SDEN,SNUM,XE,XI
      INTEGER K,M,M0,MM1
C     ..
C     .. Local Arrays ..
      REAL CDEN(8),CNUM(8)
C     ..
C     .. External Functions ..
      REAL R1MACH,T
      EXTERNAL R1MACH,T
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC AINT,INT,LOG,REAL
C     ..
C     .. Data statements ..
      DATA CNUM/4.120843185847770,85.68982062831317,243.175243524421,
     +     -261.7218583856145,-922.2613728801522,-517.6383498023218,
     +     -77.41064071332953,-2.208843997216182/,CDEN/1.,
     +     45.64677187585908,377.8372484823942,951.323597679706,
     +     846.0755362020782,262.3083470269460,24.43519662506312,
     +     .4097792921092615/
C     ..
c
c The constants in the statement below are  exp(1.)  and  .5*alog(8.).
c
      M0 = 2.71828*T((LOG(R1MACH(2))-1.03972)/2.71828)
      XI = AINT(X)
      IF (X-XI.GT..5) XI = XI + 1.
      M = INT(XI) - 1
c
c Computation of log gamma on the standard interval (1/2,3/2]
c
      XE = X - REAL(M)
      SNUM = CNUM(1)
      SDEN = CDEN(1)
      DO 10 K = 2,8
          SNUM = XE*SNUM + CNUM(K)
          SDEN = XE*SDEN + CDEN(K)
   10 CONTINUE
      ALGA = (XE-1.)*SNUM/SDEN
c
c Computation of log gamma on (0,1/2]
c
      IF (M.EQ.-1) THEN
          ALGA = ALGA - LOG(X)
          RETURN

      ELSE IF (M.EQ.0) THEN
          RETURN

      ELSE
c
c Computation of log gamma on (3/2,5/2]
c
          P = XE
          IF (M.EQ.1) THEN
              ALGA = ALGA + LOG(P)
              RETURN

          ELSE
c
c Computation of log gamma for arguments larger than 5/2
c
              MM1 = M - 1
c
c The else-clause in the next statement is designed to avoid possible
c overflow in the computation of  p  in the if-clause, at the expense
c of computing many logarithms.
c
              IF (M.LT.M0) THEN
                  DO 20 K = 1,MM1
                      P = (XE+REAL(K))*P
   20             CONTINUE
                  ALGA = ALGA + LOG(P)
                  RETURN

              ELSE
                  ALGA = ALGA + LOG(XE)
                  DO 30 K = 1,MM1
                      ALGA = ALGA + LOG(XE+REAL(K))
   30             CONTINUE
                  RETURN

              END IF

          END IF

      END IF

      END

      REAL FUNCTION GAMMA(X,IERR)
c
c This evaluates the gamma function for real positive  x, using the
c function subroutines  alga  and  r1mach. In case of overflow, the
c routine returns the largest machine-representable number and the
c error flag  ierr=2.
c
C     .. Scalar Arguments ..
      REAL X
      INTEGER IERR
C     ..
C     .. Local Scalars ..
      REAL ALMACH,T
C     ..
C     .. External Functions ..
      REAL ALGA,R1MACH
      EXTERNAL ALGA,R1MACH
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC EXP,LOG
C     ..
      ALMACH = LOG(R1MACH(2))
      IERR = 0
      T = ALGA(X)
      IF (T.GE.ALMACH) THEN
          IERR = 2
          GAMMA = R1MACH(2)
          RETURN

      ELSE
          GAMMA = EXP(T)
          RETURN

      END IF

      END

      REAL FUNCTION T(Y)
c
c This evaluates the inverse function  t = t(y)  of y = t ln t  for
c nonnegative  y  to an accuracy of about one percent. For the
c approximation used, see pp. 51-52 in W. Gautschi,Computational
c aspects of three-term recurrence relations'', SIAM Rev. 9, 1967,
c 24-82.
c
C     .. Scalar Arguments ..
      REAL Y
C     ..
C     .. Local Scalars ..
      REAL P,Z
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC LOG
C     ..
      IF (Y.LE.10.) THEN
          P = .000057941*Y - .00176148
          P = Y*P + .0208645
          P = Y*P - .129013
          P = Y*P + .85777
          T = Y*P + 1.0125

      ELSE
          Z = LOG(Y) - .775
          P = (.775-LOG(Z))/ (1.+Z)
          P = 1./ (1.+P)
          T = Y*P/Z
      END IF

      RETURN

      END
c
c
      SUBROUTINE STI(N,NCAP,X,W,ALPHA,BETA,IERR,P0,P1,P2)
c
c This routine applies Stieltjes's procedure'' (cf. Section 2.1 of
c W. Gautschi,On generating orthogonal polynomials'', SIAM J. Sci.
c Statist. Comput. 3, 1982, 289-317) to generate the recursion
c coefficients  alpha(k), beta(k) , k=0,1,...,n-1, for the discrete
c (monic) orthogonal polynomials associated with the inner product
c
c     (f,g)=sum over k from 1 to ncap of w(k)*f(x(k))*g(x(k)).
c
c The integer  n  must be between  1  and  ncap, inclusive; otherwise,
c there is an error exit with  ierr=1. The results are stored in the
c arrays  alpha, beta; the arrays  p0, p1, p2  are working arrays.
c
c If there is a threat of underflow or overflow in the calculation
c of the coefficients  alpha(k)  and  beta(k), the routine exits with
c the error flag  ierr  set equal to  -k  (in the case of underflow)
c or  +k  (in the case of overflow), where  k  is the recursion index
c for which the problem occurs. The former [latter] can often be avoided
c by multiplying all weights  w(k)  by a sufficiently large [small]
c scaling factor prior to entering the routine, and, upon exit, divide
c the coefficient  beta(0)  by the same factor.
c
c This routine should be used with caution if  n  is relatively close
c to  ncap, since there is a distinct possibility of numerical
c instability developing. (See W. Gautschi,Is the recurrence relation
c for orthogonal polynomials always stable?'', BIT, 1993, to appear.)
c In that case, the routine  lancz  should be used.
c
c The routine uses the function subroutine  r1mach.
c
C     .. Scalar Arguments ..
      INTEGER IERR,N,NCAP
C     ..
C     .. Array Arguments ..
      REAL ALPHA(N),BETA(N),P0(NCAP),P1(NCAP),P2(NCAP),W(NCAP),X(NCAP)
C     ..
C     .. Local Scalars ..
      REAL HUGE,SUM0,SUM1,SUM2,T,TINY
      INTEGER K,M,NM1
C     ..
C     .. External Functions ..
      REAL R1MACH
      EXTERNAL R1MACH
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
      TINY = 10.*R1MACH(1)
      HUGE = .1*R1MACH(2)
      IERR = 0
      IF (N.LE.0 .OR. N.GT.NCAP) THEN
          IERR = 1
          RETURN

      END IF

      NM1 = N - 1
c
c Compute the first alpha- and beta-coefficient.
c
      SUM0 = 0.
      SUM1 = 0.
      DO 10 M = 1,NCAP
          SUM0 = SUM0 + W(M)
          SUM1 = SUM1 + W(M)*X(M)
   10 CONTINUE
      ALPHA(1) = SUM1/SUM0
      BETA(1) = SUM0
      IF (N.EQ.1) RETURN
c
c Compute the remaining alpha- and beta-coefficients.
c
      DO 20 M = 1,NCAP
          P1(M) = 0.
          P2(M) = 1.
   20 CONTINUE
      DO 40 K = 1,NM1
          SUM1 = 0.
          SUM2 = 0.
          DO 30 M = 1,NCAP
c
c The following statement is designed to avoid an overflow condition
c in the computation of  p2(m)  when the weights  w(m)  go to zero
c faster (and underflow) than the  p2(m)  grow.
c
              IF (W(M).EQ.0.) GO TO 30
              P0(M) = P1(M)
              P1(M) = P2(M)
              P2(M) = (X(M)-ALPHA(K))*P1(M) - BETA(K)*P0(M)
c
c Check for impending overflow.
c
              IF (ABS(P2(M)).GT.HUGE .OR. ABS(SUM2).GT.HUGE) THEN
                  IERR = K
                  RETURN

              END IF

              T = W(M)*P2(M)*P2(M)
              SUM1 = SUM1 + T
              SUM2 = SUM2 + T*X(M)
   30     CONTINUE
c
c Check for impending underflow.
c
          IF (ABS(SUM1).LT.TINY) THEN
              IERR = -K
              RETURN

          END IF

          ALPHA(K+1) = SUM2/SUM1
          BETA(K+1) = SUM1/SUM0
          SUM0 = SUM1
   40 CONTINUE
      RETURN

      END
c
c
      INTEGER FUNCTION NU0HER(N,Z,EPS)
c
c This is an auxiliary function routine providing a starting backward
c recurrence index for the Hermite measure that can be used in place
c of  nu0  in the routines  knum  and  dknum.
c
C     .. Scalar Arguments ..
      COMPLEX Z
      REAL EPS
      INTEGER N
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,AIMAG,LOG,REAL,SQRT
C     ..
      NU0HER = 2.* (SQRT(.5*REAL(N+1))+.25*LOG(1./EPS)/ABS(AIMAG(Z)))**2
      RETURN

      END
c
c
      INTEGER FUNCTION NU0JAC(N,Z,EPS)
c
c This is an auxiliary function routine providing a starting backward
c recurrence index for the Jacobi measure that can be used in place
c of  nu0  in the routines  knum  and  dknum.
c
C     .. Scalar Arguments ..
      COMPLEX Z
      REAL EPS
      INTEGER N
C     ..
C     .. Local Scalars ..
      REAL ANGLE,PI,R,X,X2,Y,Y2
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,AIMAG,ATAN,COS,LOG,REAL,SIN,SQRT
C     ..
      PI = 4.*ATAN(1.)
      X = REAL(Z)
      Y = ABS(AIMAG(Z))
      IF (X.LT.1.) THEN
          IF (X.LT.-1.) ANGLE = .5* (2.*PI+ATAN(Y/ (X-1.))+
     +                          ATAN(Y/ (X+1.)))
          IF (X.EQ.-1.) ANGLE = .5* (1.5*PI-ATAN(.5*Y))
          IF (X.GT.-1.) ANGLE = .5* (PI+ATAN(Y/ (X-1.))+ATAN(Y/ (X+1.)))

      ELSE
          IF (X.EQ.1.) ANGLE = .5* (.5*PI+ATAN(.5*Y))
          IF (X.GT.1.) ANGLE = .5* (ATAN(Y/ (X-1.))+ATAN(Y/ (X+1.)))
      END IF

      X2 = X*X
      Y2 = Y*Y
      R = ((X2-Y2-1.)**2+4.*X2*Y2)**.25
      R = SQRT((X+R*COS(ANGLE))**2+ (Y+R*SIN(ANGLE))**2)
      NU0JAC = REAL(N+1) + .5*LOG(1./EPS)/LOG(R)
      RETURN

      END
c
c
      INTEGER FUNCTION NU0LAG(N,Z,AL,EPS)
c
c This is an auxiliary function routine providing a starting backward
c recurrence index for the Laguerre measure that can be used in place
c of  nu0  in the routines  knum  and  dknum.
c
C     .. Scalar Arguments ..
      COMPLEX Z
      REAL AL,EPS
      INTEGER N
C     ..
C     .. Local Scalars ..
      REAL PHI,PI,X,Y
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC AIMAG,ATAN,COS,LOG,REAL,SQRT
C     ..
      PI = 4.*ATAN(1.)
      X = REAL(Z)
      Y = AIMAG(Z)
      PHI = .5*PI
      IF (Y.LT.0.) PHI = 1.5*PI
      IF (X.EQ.0.) GO TO 10
      PHI = ATAN(Y/X)
      IF (Y.GT.0. .AND. X.GT.0.) GO TO 10
      PHI = PHI + PI
      IF (X.LT.0.) GO TO 10
      PHI = PHI + PI
   10 NU0LAG = (SQRT(REAL(N+1)+.5* (AL+1.))+
     +         LOG(1./EPS)/ (4.* (X*X+Y*Y)**.25*COS(.5* (PHI-PI))))**2 -
     +         .5* (AL+1.)
      RETURN

      END
