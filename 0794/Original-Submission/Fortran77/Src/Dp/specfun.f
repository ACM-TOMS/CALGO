C     BEGIN OF SUBROUTINE *** RJBESL ***
      SUBROUTINE RJBESL(XINPUT,ALPHA,NB,B,NCALC)
C
C     ##################################################################
C
C     *** T. WIEDER, 05.03.1998 ***
C
C     SOURCE FOR SUBROUTINE *** RJBESL ***:
C
C     SUBPROGRAM LIBRARY SPECFUN
C
C     AVAILABILITIY: PUBLIC DOMAIN
C
C     DEVELOPER: W.J. CODY AND L. STOLTZ, APPLIED MATHEMATICS DIVISION,
C     ARGONNE NATIONAL LABORATORY, ARGONNE, IL 60439.
C
C     DISTRIBUTOR: NETLIB
C
C     ##################################################################
C
C---------------------------------------------------------------------
C This routine calculates Bessel functions J sub(N+ALPHA) (X)
C   for non-negative argument X, and non-negative order N+ALPHA.
C
C
C  Explanation of variables in the calling sequence.
C
C   XIN = X *** T. WIEDER, 15.04.98 ***
C   X     - working precision non-negative real argument for which
C           J's are to be calculated.
C   ALPHA - working precision fractional part of order for which
C           J's or exponentially scaled J'r (J*exp(X)) are
C           to be calculated.  0 <= ALPHA < 1.0.
C   NB  - integer number of functions to be calculated, NB > 0.
C           The first function calculated is of order ALPHA, and the
C           last is of order (NB - 1 + ALPHA).
C   B  - working precision output vector of length NB.  If RJBESL
C           terminates normally (NCALC=NB), the vector B contains the
C           functions J/ALPHA/(X) through J/NB-1+ALPHA/(X), or the
C           corresponding exponentially scaled functions.
C   NCALC - integer output variable indicating possible errors.
C           Before using the vector B, the user should check that
C           NCALC=NB, i.e., all orders have been calculated to
C           the desired accuracy.  See Error Returns below.
C
C
C*******************************************************************
C*******************************************************************
C
C  Explanation of machine-dependent constants
C
C   it     = Number of bits in the mantissa of a working precision
C            variable
C   NSIG   = Decimal significance desired.  Should be set to
C            INT(LOG10(2)*it+1).  Setting NSIG lower will result
C            in decreased accuracy while setting NSIG higher will
C            increase CPU time without increasing accuracy.  The
C            truncation error is limited to a relative error of
C            T=.5*10**(-NSIG).
C   ENTEN  = 10.0 ** K, where K is the largest integer such that
C            ENTEN is machine-representable in working precision
C   ENSIG  = 10.0 ** NSIG
C   RTNSIG = 10.0 ** (-K) for the smallest integer K such that
C            K .GE. NSIG/4
C   ENMTEN = Smallest ABS(X) such that X/4 does not underflow
C   XLARGE = Upper limit on the magnitude of X.  If ABS(X)=N,
C            then at least N iterations of the backward recursion
C            will be executed.  The value of 10.0 ** 4 is used on
C            every machine.
C
C
C     Approximate values for some important machines are:
C
C
C                            it    NSIG    ENTEN       ENSIG
C
C   CRAY-1        (S.P.)     48     15    1.0E+2465   1.0E+15
C   Cyber 180/855
C     under NOS   (S.P.)     48     15    1.0E+322    1.0E+15
C   IEEE (IBM/XT,
C     SUN, etc.)  (S.P.)     24      8    1.0E+38     1.0E+8
C   IEEE (IBM/XT,
C     SUN, etc.)  (D.P.)     53     16    1.0D+308    1.0D+16
C   IBM 3033      (D.P.)     14      5    1.0D+75     1.0D+5
C   VAX           (S.P.)     24      8    1.0E+38     1.0E+8
C   VAX D-Format  (D.P.)     56     17    1.0D+38     1.0D+17
C   VAX G-Format  (D.P.)     53     16    1.0D+307    1.0D+16
C
C
C                           RTNSIG      ENMTEN      XLARGE
C
C   CRAY-1        (S.P.)    1.0E-4    1.84E-2466   1.0E+4
C   Cyber 180/855
C     under NOS   (S.P.)    1.0E-4    1.25E-293    1.0E+4
C   IEEE (IBM/XT,
C     SUN, etc.)  (S.P.)    1.0E-2    4.70E-38     1.0E+4
C   IEEE (IBM/XT,
C     SUN, etc.)  (D.P.)    1.0E-4    8.90D-308    1.0D+4
C   IBM 3033      (D.P.)    1.0E-2    2.16D-78     1.0D+4
C   VAX           (S.P.)    1.0E-2    1.17E-38     1.0E+4
C   VAX D-Format  (D.P.)    1.0E-5    1.17D-38     1.0D+4
C   VAX G-Format  (D.P.)    1.0E-4    2.22D-308    1.0D+4
C
C*******************************************************************
C*******************************************************************
C
C  Error returns
C
C    In case of an error,  NCALC .NE. NB, and not all J's are
C    calculated to the desired accuracy.
C
C    NCALC .LT. 0:  An argument is out of range. For example,
C       NBES .LE. 0, ALPHA .LT. 0 or .GT. 1, or X is too large.
C       In this case, B(1) is set to zero, the remainder of the
C       B-vector is not calculated, and NCALC is set to
C       MIN(NB,0)-1 so that NCALC .NE. NB.
C
C    NB .GT. NCALC .GT. 0: Not all requested function values could
C       be calculated accurately.  This usually occurs because NB is
C       much larger than ABS(X).  In this case, B(N) is calculated
C       to the desired accuracy for N .LE. NCALC, but precision
C       is lost for NCALC .LT. N .LE. NB.  If B(N) does not vanish
C       for N .GT. NCALC (because it is too small to be represented),
C       and B(N)/B(NCALC) = 10**(-K), then only the first NSIG-K
C       significant figures of B(N) can be trusted.
C
C
C  Intrinsic and other functions required are:
C
C     ABS, AINT, COS, DBLE, GAMMA (or DGAMMA), INT, MAX, MIN,
C
C     REAL, SIN, SQRT
C
C
C  Acknowledgement
C
C   This program is based on a program written by David J. Sookne
C   (2) that computes values of the Bessel functions J or I of real
C   argument and integer order.  Modifications include the restriction
C   of the computation to the J Bessel function of non-negative real
C   argument, the extension of the computation to arbitrary positive
C   order, and the elimination of most underflow.
C
C  References: "A Note on Backward Recurrence Algorithms," Olver,
C               F. W. J., and Sookne, D. J., Math. Comp. 26, 1972,
C               pp 941-947.
C
C              "Bessel Functions of Real Argument and Integer Order,"
C               Sookne, D. J., NBS Jour. of Res. B. 77B, 1973, pp
C               125-132.
C
C  Latest modification: March 19, 1990
C
C  Author: W. J. Cody
C          Applied Mathematics Division
C          Argonne National Laboratory
C          Argonne, IL  60439
C
C---------------------------------------------------------------------
CS    REAL               GAMMA,
C     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA,XINPUT
      INTEGER NB,NCALC
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION B(NB)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ALP2EM,ALPEM,CAPP,CAPQ,EIGHTH,EM,EN,ENMTEN,ENSIG,
     +                 ENTEN,FOUR,GNU,HALF,HALFX,ONE,ONE30,P,PI2,PLAST,
     +                 POLD,PSAVE,PSAVEL,RTNSIG,S,SUM,T,T1,TEMPA,TEMPB,
     +                 TEMPC,TEST,THREE,THREE5,TOVER,TWO,TWOFIV,TWOPI1,
     +                 TWOPI2,VCOS,VSIN,X,XC,XIN,XK,XLARGE,XM,Z,ZERO
      INTEGER I,J,K,L,M,MAGX,N,NBMX,NEND,NSTART
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION FACT(25)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DGAMMA
      EXTERNAL DGAMMA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,AINT,COS,DBLE,INT,MAX,MIN,SIN,SQRT
C     ..
C     .. Statement Functions ..
      DOUBLE PRECISION CONV,FUNC
C     ..
C     .. Data statements ..
C---------------------------------------------------------------------
C  Mathematical constants
C
C   PI2    - 2 / PI
C   TWOPI1 - first few significant digits of 2 * PI
C   TWOPI2 - (2*PI - TWOPI) to working precision, i.e.,
C            TWOPI1 + TWOPI2 = 2 * PI to extra precision.
C---------------------------------------------------------------------
C    DATA PI2, TWOPI1, TWOPI2 /0.636619772367581343075535E0,6.28125E0,
C   1 1.935307179586476925286767E-3/
C    DATA ZERO, EIGHTH, HALF, ONE /0.0E0,0.125E0,0.5E0,1.0E0/
C    DATA TWO, THREE, FOUR, TWOFIV /2.0E0,3.0E0,4.0E0,25.0E0/
C    DATA ONE30, THREE5 /130.0E0,35.0E0/
C---------------------------------------------------------------------
C  Machine-dependent parameters
C---------------------------------------------------------------------
CS    DATA ENTEN, ENSIG, RTNSIG /1.0E38,1.0E8,1.0E-2/
CS    DATA ENMTEN, XLARGE /1.2E-37,1.0E4/
C---------------------------------------------------------------------
C     Factorial(N)
C---------------------------------------------------------------------
CS    DATA FACT /1.0E0,1.0E0,2.0E0,6.0E0,24.0E0,1.2E2,7.2E2,5.04E3,
CS   1 4.032E4,3.6288E5,3.6288E6,3.99168E7,4.790016E8,6.2270208E9,
CS   2 8.71782912E10,1.307674368E12,2.0922789888E13,3.55687428096E14,
CS   3 6.402373705728E15,1.21645100408832E17,2.43290200817664E18,
CS   4 5.109094217170944E19,1.12400072777760768E21,
CS   5 2.585201673888497664E22,6.2044840173323943936E23/
      DATA PI2,TWOPI1,TWOPI2/0.636619772367581343075535D0,6.28125D0,
     +     1.935307179586476925286767D-3/
      DATA ZERO,EIGHTH,HALF,ONE/0.0D0,0.125D0,0.5D0,1.0D0/
      DATA TWO,THREE,FOUR,TWOFIV/2.0D0,3.0D0,4.0D0,25.0D0/
      DATA ONE30,THREE5/130.0D0,35.0D0/
      DATA ENTEN,ENSIG,RTNSIG/1.0D38,1.0D17,1.0D-4/
      DATA ENMTEN,XLARGE/1.2D-37,1.0D4/
      DATA FACT/1.0D0,1.0D0,2.0D0,6.0D0,24.0D0,1.2D2,7.2D2,5.04D3,
     +     4.032D4,3.6288D5,3.6288D6,3.99168D7,4.790016D8,6.2270208D9,
     +     8.71782912D10,1.307674368D12,2.0922789888D13,
     +     3.55687428096D14,6.402373705728D15,1.21645100408832D17,
     +     2.43290200817664D18,5.109094217170944D19,
     +     1.12400072777760768D21,2.585201673888497664D22,
     +     6.2044840173323943936D23/
C     ..
C     .. Statement Function definitions ..
C
C---------------------------------------------------------------------
C Statement functions for conversion and the gamma function.
C---------------------------------------------------------------------
CS    CONV(I) = REAL(I)
CS    FUNC(X) = GAMMA(X)
      CONV(I) = DBLE(I)
      FUNC(X) = DGAMMA(X)
C     ..
C
C     *** T. WIEDER, 15.04.1998 ***
C     *** AMENDMENT SUGGESTED BY DR. T.R. HOPKINS, 15.04.1998 ***
      X = XINPUT
C
C---------------------------------------------------------------------
C Check for out of range arguments.
C---------------------------------------------------------------------
      MAGX = INT(X)
      IF ((NB.GT.0) .AND. (X.GE.ZERO) .AND. (X.LE.XLARGE) .AND.
     +    (ALPHA.GE.ZERO) .AND. (ALPHA.LT.ONE)) THEN
C---------------------------------------------------------------------
C Initialize result array to zero.
C---------------------------------------------------------------------
          NCALC = NB
          DO 10 I = 1,NB
              B(I) = ZERO
   10     CONTINUE
C---------------------------------------------------------------------
C Branch to use 2-term ascending series for small X and asymptotic
C form for large X when NB is not too large.
C---------------------------------------------------------------------
          IF (X.LT.RTNSIG) THEN
C---------------------------------------------------------------------
C Two-term ascending series for small X.
C---------------------------------------------------------------------
              TEMPA = ONE
              ALPEM = ONE + ALPHA
              HALFX = ZERO
              IF (X.GT.ENMTEN) HALFX = HALF*X
              IF (ALPHA.NE.ZERO) TEMPA = HALFX**ALPHA/
     +                                   (ALPHA*FUNC(ALPHA))
              TEMPB = ZERO
              IF ((X+ONE).GT.ONE) TEMPB = -HALFX*HALFX
              B(1) = TEMPA + TEMPA*TEMPB/ALPEM
              IF ((X.NE.ZERO) .AND. (B(1).EQ.ZERO)) NCALC = 0
              IF (NB.NE.1) THEN
                  IF (X.LE.ZERO) THEN
                      DO 20 N = 2,NB
                          B(N) = ZERO
   20                 CONTINUE

                  ELSE
C---------------------------------------------------------------------
C Calculate higher order functions.
C---------------------------------------------------------------------
                      TEMPC = HALFX
                      TOVER = (ENMTEN+ENMTEN)/X
                      IF (TEMPB.NE.ZERO) TOVER = ENMTEN/TEMPB
                      DO 30 N = 2,NB
                          TEMPA = TEMPA/ALPEM
                          ALPEM = ALPEM + ONE
                          TEMPA = TEMPA*TEMPC
                          IF (TEMPA.LE.TOVER*ALPEM) TEMPA = ZERO
                          B(N) = TEMPA + TEMPA*TEMPB/ALPEM
                          IF ((B(N).EQ.ZERO) .AND.
     +                        (NCALC.GT.N)) NCALC = N - 1
   30                 CONTINUE
                  END IF

              END IF

          ELSE IF ((X.GT.TWOFIV) .AND. (NB.LE.MAGX+1)) THEN
C---------------------------------------------------------------------
C Asymptotic series for X .GT. 21.0.
C---------------------------------------------------------------------
              XC = SQRT(PI2/X)
              XIN = (EIGHTH/X)**2
              M = 11
              IF (X.GE.THREE5) M = 8
              IF (X.GE.ONE30) M = 4
              XM = FOUR*CONV(M)
C---------------------------------------------------------------------
C Argument reduction for SIN and COS routines.
C---------------------------------------------------------------------
              T = AINT(X/ (TWOPI1+TWOPI2)+HALF)
              Z = ((X-T*TWOPI1)-T*TWOPI2) - (ALPHA+HALF)/PI2
              VSIN = SIN(Z)
              VCOS = COS(Z)
              GNU = ALPHA + ALPHA
              DO 50 I = 1,2
                  S = ((XM-ONE)-GNU)* ((XM-ONE)+GNU)*XIN*HALF
                  T = (GNU- (XM-THREE))* (GNU+ (XM-THREE))
                  CAPP = S*T/FACT(2*M+1)
                  T1 = (GNU- (XM+ONE))* (GNU+ (XM+ONE))
                  CAPQ = S*T1/FACT(2*M+2)
                  XK = XM
                  K = M + M
                  T1 = T
                  DO 40 J = 2,M
                      XK = XK - FOUR
                      S = ((XK-ONE)-GNU)* ((XK-ONE)+GNU)
                      T = (GNU- (XK-THREE))* (GNU+ (XK-THREE))
                      CAPP = (CAPP+ONE/FACT(K-1))*S*T*XIN
                      CAPQ = (CAPQ+ONE/FACT(K))*S*T1*XIN
                      K = K - 2
                      T1 = T
   40             CONTINUE
                  CAPP = CAPP + ONE
                  CAPQ = (CAPQ+ONE)* (GNU*GNU-ONE)* (EIGHTH/X)
                  B(I) = XC* (CAPP*VCOS-CAPQ*VSIN)
                  IF (NB.EQ.1) GO TO 180
                  T = VSIN
                  VSIN = -VCOS
                  VCOS = T
                  GNU = GNU + TWO
   50         CONTINUE
C---------------------------------------------------------------------
C If  NB .GT. 2, compute J(X,ORDER+I)  I = 2, NB-1
C---------------------------------------------------------------------
              IF (NB.GT.2) THEN
                  GNU = ALPHA + ALPHA + TWO
                  DO 60 J = 3,NB
                      B(J) = GNU*B(J-1)/X - B(J-2)
                      GNU = GNU + TWO
   60             CONTINUE
              END IF
C---------------------------------------------------------------------
C Use recurrence to generate results.  First initialize the
C calculation of P*S.
C---------------------------------------------------------------------
          ELSE
              NBMX = NB - MAGX
              N = MAGX + 1
              EN = CONV(N+N) + (ALPHA+ALPHA)
              PLAST = ONE
              P = EN/X
C---------------------------------------------------------------------
C Calculate general significance test.
C---------------------------------------------------------------------
              TEST = ENSIG + ENSIG
              IF (NBMX.GE.3) THEN
C---------------------------------------------------------------------
C Calculate P*S until N = NB-1.  Check for possible overflow.
C---------------------------------------------------------------------
                  TOVER = ENTEN/ENSIG
                  NSTART = MAGX + 2
                  NEND = NB - 1
                  EN = CONV(NSTART+NSTART) - TWO + (ALPHA+ALPHA)
                  DO 90 K = NSTART,NEND
                      N = K
                      EN = EN + TWO
                      POLD = PLAST
                      PLAST = P
                      P = EN*PLAST/X - POLD
                      IF (P.GT.TOVER) THEN
C---------------------------------------------------------------------
C To avoid overflow, divide P*S by TOVER.  Calculate P*S until
C ABS(P) .GT. 1.
C---------------------------------------------------------------------
                          TOVER = ENTEN
                          P = P/TOVER
                          PLAST = PLAST/TOVER
                          PSAVE = P
                          PSAVEL = PLAST
                          NSTART = N + 1
   70                     N = N + 1
                          EN = EN + TWO
                          POLD = PLAST
                          PLAST = P
                          P = EN*PLAST/X - POLD
                          IF (P.LE.ONE) GO TO 70
                          TEMPB = EN/X
C---------------------------------------------------------------------
C Calculate backward test and find NCALC, the highest N such that
C the test is passed.
C---------------------------------------------------------------------
                          TEST = POLD*PLAST* (HALF-HALF/ (TEMPB*TEMPB))
                          TEST = TEST/ENSIG
                          P = PLAST*TOVER
                          N = N - 1
                          EN = EN - TWO
                          NEND = MIN(NB,N)
                          DO 80 L = NSTART,NEND
                              POLD = PSAVEL
                              PSAVEL = PSAVE
                              PSAVE = EN*PSAVEL/X - POLD
                              IF (PSAVE*PSAVEL.GT.TEST) THEN
                                  NCALC = L - 1
                                  GO TO 110

                              END IF

   80                     CONTINUE
                          NCALC = NEND
                          GO TO 110

                      END IF

   90             CONTINUE
                  N = NEND
                  EN = CONV(N+N) + (ALPHA+ALPHA)
C---------------------------------------------------------------------
C Calculate special significance test for NBMX .GT. 2.
C---------------------------------------------------------------------
                  TEST = MAX(TEST,SQRT(PLAST*ENSIG)*SQRT(P+P))
              END IF
C---------------------------------------------------------------------
C Calculate P*S until significance test passes.
C---------------------------------------------------------------------
  100         N = N + 1
              EN = EN + TWO
              POLD = PLAST
              PLAST = P
              P = EN*PLAST/X - POLD
              IF (P.LT.TEST) GO TO 100
C---------------------------------------------------------------------
C Initialize the backward recursion and the normalization sum.
C---------------------------------------------------------------------
  110         N = N + 1
              EN = EN + TWO
              TEMPB = ZERO
              TEMPA = ONE/P
              M = 2*N - 4* (N/2)
              SUM = ZERO
              EM = CONV(N/2)
              ALPEM = (EM-ONE) + ALPHA
              ALP2EM = (EM+EM) + ALPHA
              IF (M.NE.0) SUM = TEMPA*ALPEM*ALP2EM/EM
              NEND = N - NB
              IF (NEND.GT.0) THEN
C---------------------------------------------------------------------
C Recur backward via difference equation, calculating (but not
C storing) B(N), until N = NB.
C---------------------------------------------------------------------
                  DO 120 L = 1,NEND
                      N = N - 1
                      EN = EN - TWO
                      TEMPC = TEMPB
                      TEMPB = TEMPA
                      TEMPA = (EN*TEMPB)/X - TEMPC
                      M = 2 - M
                      IF (M.NE.0) THEN
                          EM = EM - ONE
                          ALP2EM = (EM+EM) + ALPHA
                          IF (N.EQ.1) GO TO 130
                          ALPEM = (EM-ONE) + ALPHA
                          IF (ALPEM.EQ.ZERO) ALPEM = ONE
                          SUM = (SUM+TEMPA*ALP2EM)*ALPEM/EM
                      END IF

  120             CONTINUE
              END IF
C---------------------------------------------------------------------
C Store B(NB).
C---------------------------------------------------------------------
  130         B(N) = TEMPA
              IF (NEND.GE.0) THEN
                  IF (NB.LE.1) THEN
                      ALP2EM = ALPHA
                      IF ((ALPHA+ONE).EQ.ONE) ALP2EM = ONE
                      SUM = SUM + B(1)*ALP2EM
                      GO TO 160

                  ELSE
C---------------------------------------------------------------------
C Calculate and store B(NB-1).
C---------------------------------------------------------------------
                      N = N - 1
                      EN = EN - TWO
                      B(N) = (EN*TEMPA)/X - TEMPB
                      IF (N.EQ.1) GO TO 150
                      M = 2 - M
                      IF (M.NE.0) THEN
                          EM = EM - ONE
                          ALP2EM = (EM+EM) + ALPHA
                          ALPEM = (EM-ONE) + ALPHA
                          IF (ALPEM.EQ.ZERO) ALPEM = ONE
                          SUM = (SUM+B(N)*ALP2EM)*ALPEM/EM
                      END IF

                  END IF

              END IF

              NEND = N - 2
              IF (NEND.NE.0) THEN
C---------------------------------------------------------------------
C Calculate via difference equation and store B(N), until N = 2.
C---------------------------------------------------------------------
                  DO 140 L = 1,NEND
                      N = N - 1
                      EN = EN - TWO
                      B(N) = (EN*B(N+1))/X - B(N+2)
                      M = 2 - M
                      IF (M.NE.0) THEN
                          EM = EM - ONE
                          ALP2EM = (EM+EM) + ALPHA
                          ALPEM = (EM-ONE) + ALPHA
                          IF (ALPEM.EQ.ZERO) ALPEM = ONE
                          SUM = (SUM+B(N)*ALP2EM)*ALPEM/EM
                      END IF

  140             CONTINUE
              END IF
C---------------------------------------------------------------------
C Calculate B(1).
C---------------------------------------------------------------------
              B(1) = TWO* (ALPHA+ONE)*B(2)/X - B(3)
  150         EM = EM - ONE
              ALP2EM = (EM+EM) + ALPHA
              IF (ALP2EM.EQ.ZERO) ALP2EM = ONE
              SUM = SUM + B(1)*ALP2EM
C---------------------------------------------------------------------
C Normalize.  Divide all B(N) by sum.
C---------------------------------------------------------------------
  160         IF ((ALPHA+ONE).NE.ONE) SUM = SUM*FUNC(ALPHA)*
     +            (X*HALF)** (-ALPHA)
              TEMPA = ENMTEN
              IF (SUM.GT.ONE) TEMPA = TEMPA*SUM
              DO 170 N = 1,NB
                  IF (ABS(B(N)).LT.TEMPA) B(N) = ZERO
                  B(N) = B(N)/SUM
  170         CONTINUE
          END IF
C---------------------------------------------------------------------
C Error return -- X, NB, or ALPHA is out of range.
C---------------------------------------------------------------------
      ELSE
          B(1) = ZERO
          NCALC = MIN(NB,0) - 1
      END IF
C---------------------------------------------------------------------
C Exit
C---------------------------------------------------------------------
  180 RETURN
C ---------- Last line of RJBESL ----------
      END
C     END OF SUBROUTINE *** RJBESL ***
C
C     BEGIN OF SUBROUTINE *** RYBESL ***
      SUBROUTINE RYBESL(X,ALPHA,NB,BY,NCALC)
C
C     ##################################################################
C
C     *** T. WIEDER, 05.03.1998 ***
C
C     SOURCE FOR *** SUBROUTINE RYBESL ***:
C
C     SUBPROGRAM LIBRARY SPECFUN
C
C     AVAILABILITIY: PUBLIC DOMAIN
C
C     DEVELOPER: W.J. CODY AND L. STOLTZ, APPLIED MATHEMATICS DIVISION,
C     ARGONNE NATIONAL LABORATORY, ARGONNE, IL 60439.
C
C     DISTRIBUTOR: NETLIB
C
C     ##################################################################
C
C----------------------------------------------------------------------
C
C  This routine calculates Bessel functions Y SUB(N+ALPHA) (X)
C  for non-negative argument X, and non-negative order N+ALPHA.
C
C
C Explanation of variables in the calling sequence
C
C X     - Working precision non-negative real argument for which
C         Y's are to be calculated.
C ALPHA - Working precision fractional part of order for which
C         Y's are to be calculated.  0 .LE. ALPHA .LT. 1.0.
C NB    - Integer number of functions to be calculated, NB .GT. 0.
C         The first function calculated is of order ALPHA, and the
C         last is of order (NB - 1 + ALPHA).
C BY    - Working precision output vector of length NB.  If the
C         routine terminates normally (NCALC=NB), the vector BY
C         contains the functions Y(ALPHA,X), ... , Y(NB-1+ALPHA,X),
C         If (0 .LT. NCALC .LT. NB), BY(I) contains correct function
C         values for I .LE. NCALC, and contains the ratios
C         Y(ALPHA+I-1,X)/Y(ALPHA+I-2,X) for the rest of the array.
C NCALC - Integer output variable indicating possible errors.
C         Before using the vector BY, the user should check that
C         NCALC=NB, i.e., all orders have been calculated to
C         the desired accuracy.  See error returns below.
C
C
C*******************************************************************
C*******************************************************************
C
C Explanation of machine-dependent constants
C
C   beta   = Radix for the floating-point system
C   p      = Number of significant base-beta digits in the
C            significand of a floating-point number
C   minexp = Smallest representable power of beta
C   maxexp = Smallest power of beta that overflows
C   EPS    = beta ** (-p)
C   DEL    = Machine number below which sin(x)/x = 1; approximately
C            SQRT(EPS).
C   XMIN   = Smallest acceptable argument for RBESY; approximately
C            max(2*beta**minexp,2/XINF), rounded up
C   XINF   = Largest positive machine number; approximately
C            beta**maxexp
C   THRESH = Lower bound for use of the asymptotic form; approximately
C            AINT(-LOG10(EPS/2.0))+1.0
C   XLARGE = Upper bound on X; approximately 1/DEL, because the sine
C            and cosine functions have lost about half of their
C            precision at that point.
C
C
C     Approximate values for some important machines are:
C
C                        beta    p     minexp      maxexp      EPS
C
C  CRAY-1        (S.P.)    2    48     -8193        8191    3.55E-15
C  Cyber 180/185
C    under NOS   (S.P.)    2    48      -975        1070    3.55E-15
C  IEEE (IBM/XT,
C    SUN, etc.)  (S.P.)    2    24      -126         128    5.96E-8
C  IEEE (IBM/XT,
C    SUN, etc.)  (D.P.)    2    53     -1022        1024    1.11D-16
C  IBM 3033      (D.P.)   16    14       -65          63    1.39D-17
C  VAX           (S.P.)    2    24      -128         127    5.96E-8
C  VAX D-Format  (D.P.)    2    56      -128         127    1.39D-17
C  VAX G-Format  (D.P.)    2    53     -1024        1023    1.11D-16
C
C
C                         DEL      XMIN      XINF     THRESH  XLARGE
C
C CRAY-1        (S.P.)  5.0E-8  3.67E-2466 5.45E+2465  15.0E0  2.0E7
C Cyber 180/855
C   under NOS   (S.P.)  5.0E-8  6.28E-294  1.26E+322   15.0E0  2.0E7
C IEEE (IBM/XT,
C   SUN, etc.)  (S.P.)  1.0E-4  2.36E-38   3.40E+38     8.0E0  1.0E4
C IEEE (IBM/XT,
C   SUN, etc.)  (D.P.)  1.0D-8  4.46D-308  1.79D+308   16.0D0  1.0D8
C IBM 3033      (D.P.)  1.0D-8  2.77D-76   7.23D+75    17.0D0  1.0D8
C VAX           (S.P.)  1.0E-4  1.18E-38   1.70E+38     8.0E0  1.0E4
C VAX D-Format  (D.P.)  1.0D-9  1.18D-38   1.70D+38    17.0D0  1.0D9
C VAX G-Format  (D.P.)  1.0D-8  2.23D-308  8.98D+307   16.0D0  1.0D8
C
C*******************************************************************
C*******************************************************************
C
C Error returns
C
C  In case of an error, NCALC .NE. NB, and not all Y's are
C  calculated to the desired accuracy.
C
C  NCALC .LT. -1:  An argument is out of range. For example,
C       NB .LE. 0, IZE is not 1 or 2, or IZE=1 and ABS(X) .GE.
C       XMAX.  In this case, BY(1) = 0.0, the remainder of the
C       BY-vector is not calculated, and NCALC is set to
C       MIN0(NB,0)-2  so that NCALC .NE. NB.
C  NCALC = -1:  Y(ALPHA,X) .GE. XINF.  The requested function
C       values are set to 0.0.
C  1 .LT. NCALC .LT. NB: Not all requested function values could
C       be calculated accurately.  BY(I) contains correct function
C       values for I .LE. NCALC, and and the remaining NB-NCALC
C       array elements contain 0.0.
C
C
C Intrinsic functions required are:
C
C     DBLE, EXP, INT, MAX, MIN, REAL, SQRT
C
C
C Acknowledgement
C
C  This program draws heavily on Temme's Algol program for Y(a,x)
C  and Y(a+1,x) and on Campbell's programs for Y_nu(x).  Temme's
C  scheme is used for  x < THRESH, and Campbell's scheme is used
C  in the asymptotic region.  Segments of code from both sources
C  have been translated into Fortran 77, merged, and heavily modified.
C  Modifications include parameterization of machine dependencies,
C  use of a new approximation for ln(gamma(x)), and built-in
C  protection against over/underflow.
C
C References: "Bessel functions J_nu(x) and Y_nu(x) of real
C              order and real argument," Campbell, J. B.,
C              Comp. Phy. Comm. 18, 1979, pp. 133-142.
C
C             "On the numerical evaluation of the ordinary
C              Bessel function of the second kind," Temme,
C              N. M., J. Comput. Phys. 21, 1976, pp. 343-350.
C
C  Latest modification: March 19, 1990
C
C  Modified by: W. J. Cody
C               Applied Mathematics Division
C               Argonne National Laboratory
C               Argonne, IL  60439
C
C----------------------------------------------------------------------
CS    REAL
C      DIMENSION BY(NB), CH(21)
C     *** T. WIEDER, 07.02.1999 ***
C     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA,X
      INTEGER NB,NCALC
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION BY(NB+1)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ALFA,AYE,B,C,COSMU,D,D1,D2,DDIV,DEL,DEN,DIV,DMU,
     +                 E,EIGHT,EN,EN1,ENU,EPS,EVEN,EX,F,FIVPI,G,GAMMA,H,
     +                 HALF,ODD,ONBPI,ONE,ONE5,P,PA,PA1,PI,PIBY2,PIM5,Q,
     +                 Q0,QA,QA1,R,S,SINMU,SQ2BPI,TEN9,TERM,THREE,
     +                 THRESH,TWO,TWOBYX,X2,XINF,XLARGE,XMIN,XNA,YA,YA1,
     +                 ZERO
      INTEGER I,K,NA
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION CH(21)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,AINT,COS,INT,LOG,MIN,SIN,SQRT
C     ..
C     .. Data statements ..
C----------------------------------------------------------------------
C  Mathematical constants
C    FIVPI = 5*PI
C    PIM5 = 5*PI - 15
C    ONBPI = 1/PI
C    PIBY2 = PI/2
C    SQ2BPI = SQUARE ROOT OF 2/PI
C----------------------------------------------------------------------
CS    DATA ZERO,HALF,ONE,TWO,THREE/0.0E0,0.5E0,1.0E0,2.0E0,3.0E0/
CS    DATA EIGHT,ONE5,TEN9/8.0E0,15.0E0,1.9E1/
CS    DATA FIVPI,PIBY2/1.5707963267948966192E1,1.5707963267948966192E0/
CS    DATA PI,SQ2BPI/3.1415926535897932385E0,7.9788456080286535588E-1/
CS    DATA PIM5,ONBPI/7.0796326794896619231E-1,3.1830988618379067154E-1/
C----------------------------------------------------------------------
C  Machine-dependent constants
C----------------------------------------------------------------------
CS    DATA DEL,XMIN,XINF,EPS/1.0E-4,2.36E-38,3.40E38,5.96E-8/
CS    DATA THRESH,XLARGE/8.0E0,1.0E4/
C     *** T. WIEDER, 11.03.1998 ***
C      DATA DEL,XMIN,XINF,EPS/1.0D-8,4.46D-308,1.79D308,1.11D-16/
C----------------------------------------------------------------------
C  Coefficients for Chebyshev polynomial expansion of
C         1/gamma(1-x), abs(x) .le. .5
C----------------------------------------------------------------------
CS    DATA CH/-0.67735241822398840964E-23,-0.61455180116049879894E-22,
CS   1         0.29017595056104745456E-20, 0.13639417919073099464E-18,
CS   2         0.23826220476859635824E-17,-0.90642907957550702534E-17,
CS   3        -0.14943667065169001769E-14,-0.33919078305362211264E-13,
CS   4        -0.17023776642512729175E-12, 0.91609750938768647911E-11,
CS   5         0.24230957900482704055E-09, 0.17451364971382984243E-08,
CS   6        -0.33126119768180852711E-07,-0.86592079961391259661E-06,
CS   7        -0.49717367041957398581E-05, 0.76309597585908126618E-04,
CS   8         0.12719271366545622927E-02, 0.17063050710955562222E-02,
CS   9        -0.76852840844786673690E-01,-0.28387654227602353814E+00,
CS   A         0.92187029365045265648E+00/
      DATA ZERO,HALF,ONE,TWO,THREE/0.0D0,0.5D0,1.0D0,2.0D0,3.0D0/
      DATA EIGHT,ONE5,TEN9/8.0D0,15.0D0,1.9D1/
      DATA FIVPI,PIBY2/1.5707963267948966192D1,1.5707963267948966192D0/
      DATA PI,SQ2BPI/3.1415926535897932385D0,7.9788456080286535588D-1/
      DATA PIM5,ONBPI/7.0796326794896619231D-1,3.1830988618379067154D-1/
      DATA DEL,XMIN,XINF,EPS/1.0D-8,2.23D-307,1.79D308,1.11D-16/
      DATA THRESH,XLARGE/16.0D0,1.0D8/
      DATA CH/-0.67735241822398840964D-23,-0.61455180116049879894D-22,
     +     0.29017595056104745456D-20,0.13639417919073099464D-18,
     +     0.23826220476859635824D-17,-0.90642907957550702534D-17,
     +     -0.14943667065169001769D-14,-0.33919078305362211264D-13,
     +     -0.17023776642512729175D-12,0.91609750938768647911D-11,
     +     0.24230957900482704055D-09,0.17451364971382984243D-08,
     +     -0.33126119768180852711D-07,-0.86592079961391259661D-06,
     +     -0.49717367041957398581D-05,0.76309597585908126618D-04,
     +     0.12719271366545622927D-02,0.17063050710955562222D-02,
     +     -0.76852840844786673690D-01,-0.28387654227602353814D+00,
     +     0.92187029365045265648D+00/
C     ..
C----------------------------------------------------------------------
      EX = X
      ENU = ALPHA
      IF ((NB.GT.0) .AND. (X.GE.XMIN) .AND. (EX.LT.XLARGE) .AND.
     +    (ENU.GE.ZERO) .AND. (ENU.LT.ONE)) THEN
          XNA = AINT(ENU+HALF)
          NA = INT(XNA)
          IF (NA.EQ.1) ENU = ENU - XNA
          IF (ENU.EQ.-HALF) THEN
              P = SQ2BPI/SQRT(EX)
              YA = P*SIN(EX)
              YA1 = -P*COS(EX)

          ELSE IF (EX.LT.THREE) THEN
C----------------------------------------------------------------------
C  Use Temme's scheme for small X
C----------------------------------------------------------------------
              B = EX*HALF
              D = -LOG(B)
              F = ENU*D
              E = B** (-ENU)
              IF (ABS(ENU).LT.DEL) THEN
                  C = ONBPI

              ELSE
                  C = ENU/SIN(ENU*PI)
              END IF
C----------------------------------------------------------------------
C  Computation of sinh(f)/f
C----------------------------------------------------------------------
              IF (ABS(F).LT.ONE) THEN
                  X2 = F*F
                  EN = TEN9
                  S = ONE
                  DO 10 I = 1,9
                      S = S*X2/EN/ (EN-ONE) + ONE
                      EN = EN - TWO
   10             CONTINUE

              ELSE
                  S = (E-ONE/E)*HALF/F
              END IF
C----------------------------------------------------------------------
C  Computation of 1/gamma(1-a) using Chebyshev polynomials
C----------------------------------------------------------------------
              X2 = ENU*ENU*EIGHT
              AYE = CH(1)
              EVEN = ZERO
              ALFA = CH(2)
              ODD = ZERO
              DO 20 I = 3,19,2
                  EVEN = - (AYE+AYE+EVEN)
                  AYE = -EVEN*X2 - AYE + CH(I)
                  ODD = - (ALFA+ALFA+ODD)
                  ALFA = -ODD*X2 - ALFA + CH(I+1)
   20         CONTINUE
              EVEN = (EVEN*HALF+AYE)*X2 - AYE + CH(21)
              ODD = (ODD+ALFA)*TWO
              GAMMA = ODD*ENU + EVEN
C----------------------------------------------------------------------
C  End of computation of 1/gamma(1-a)
C----------------------------------------------------------------------
              G = E*GAMMA
              E = (E+ONE/E)*HALF
              F = TWO*C* (ODD*E+EVEN*S*D)
              E = ENU*ENU
              P = G*C
              Q = ONBPI/G
              C = ENU*PIBY2
              IF (ABS(C).LT.DEL) THEN
                  R = ONE

              ELSE
                  R = SIN(C)/C
              END IF

              R = PI*C*R*R
              C = ONE
              D = -B*B
              H = ZERO
              YA = F + R*Q
              YA1 = P
              EN = ZERO
   30         EN = EN + ONE
              IF (ABS(G/ (ONE+ABS(YA)))+ABS(H/ (ONE+ABS(YA1))).GT.
     +            EPS) THEN
                  F = (F*EN+P+Q)/ (EN*EN-E)
                  C = C*D/EN
                  P = P/ (EN-ENU)
                  Q = Q/ (EN+ENU)
                  G = C* (F+R*Q)
                  H = C*P - EN*G
                  YA = YA + G
                  YA1 = YA1 + H
                  GO TO 30

              END IF

              YA = -YA
              YA1 = -YA1/B

          ELSE IF (EX.LT.THRESH) THEN
C----------------------------------------------------------------------
C  Use Temme's scheme for moderate X
C----------------------------------------------------------------------
              C = (HALF-ENU)* (HALF+ENU)
              B = EX + EX
              E = (EX*ONBPI*COS(ENU*PI)/EPS)
              E = E*E
              P = ONE
              Q = -EX
              R = ONE + EX*EX
              S = R
              EN = TWO
   40         IF (R*EN*EN.LT.E) THEN
                  EN1 = EN + ONE
                  D = (EN-ONE+C/EN)/S
                  P = (EN+EN-P*D)/EN1
                  Q = (-B+Q*D)/EN1
                  S = P*P + Q*Q
                  R = R*S
                  EN = EN1
                  GO TO 40

              END IF

              F = P/S
              P = F
              G = -Q/S
              Q = G
   50         EN = EN - ONE
              IF (EN.GT.ZERO) THEN
                  R = EN1* (TWO-P) - TWO
                  S = B + EN1*Q
                  D = (EN-ONE+C/EN)/ (R*R+S*S)
                  P = D*R
                  Q = D*S
                  E = F + ONE
                  F = P*E - G*Q
                  G = Q*E + P*G
                  EN1 = EN
                  GO TO 50

              END IF

              F = ONE + F
              D = F*F + G*G
              PA = F/D
              QA = -G/D
              D = ENU + HALF - P
              Q = Q + EX
              PA1 = (PA*Q-QA*D)/EX
              QA1 = (QA*Q+PA*D)/EX
              B = EX - PIBY2* (ENU+HALF)
              C = COS(B)
              S = SIN(B)
              D = SQ2BPI/SQRT(EX)
              YA = D* (PA*S+QA*C)
              YA1 = D* (QA1*S-PA1*C)

          ELSE
C----------------------------------------------------------------------
C  Use Campbell's asymptotic scheme.
C----------------------------------------------------------------------
              NA = 0
              D1 = AINT(EX/FIVPI)
              I = INT(D1)
              DMU = ((EX-ONE5*D1)-D1*PIM5) - (ALPHA+HALF)*PIBY2
              IF (I-2* (I/2).EQ.0) THEN
                  COSMU = COS(DMU)
                  SINMU = SIN(DMU)

              ELSE
                  COSMU = -COS(DMU)
                  SINMU = -SIN(DMU)
              END IF

              DDIV = EIGHT*EX
              DMU = ALPHA
              DEN = SQRT(EX)
              DO 80 K = 1,2
                  P = COSMU
                  COSMU = SINMU
                  SINMU = -P
                  D1 = (TWO*DMU-ONE)* (TWO*DMU+ONE)
                  D2 = ZERO
                  DIV = DDIV
                  P = ZERO
                  Q = ZERO
                  Q0 = D1/DIV
                  TERM = Q0
                  DO 60 I = 2,20
                      D2 = D2 + EIGHT
                      D1 = D1 - D2
                      DIV = DIV + DDIV
                      TERM = -TERM*D1/DIV
                      P = P + TERM
                      D2 = D2 + EIGHT
                      D1 = D1 - D2
                      DIV = DIV + DDIV
                      TERM = TERM*D1/DIV
                      Q = Q + TERM
                      IF (ABS(TERM).LE.EPS) GO TO 70
   60             CONTINUE
   70             P = P + ONE
                  Q = Q + Q0
                  IF (K.EQ.1) THEN
                      YA = SQ2BPI* (P*COSMU-Q*SINMU)/DEN

                  ELSE
                      YA1 = SQ2BPI* (P*COSMU-Q*SINMU)/DEN
                  END IF

                  DMU = DMU + ONE
   80         CONTINUE
          END IF

          IF (NA.EQ.1) THEN
              H = TWO* (ENU+ONE)/EX
              IF (H.GT.ONE) THEN
                  IF (ABS(YA1).GT.XINF/H) THEN
                      H = ZERO
                      YA = ZERO
                  END IF

              END IF

              H = H*YA1 - YA
              YA = YA1
              YA1 = H
          END IF
C----------------------------------------------------------------------
C  Now have first one or two Y's
C----------------------------------------------------------------------
          BY(1) = YA
          BY(2) = YA1
          IF (YA1.EQ.ZERO) THEN
              NCALC = 1

          ELSE
              AYE = ONE + ALPHA
              TWOBYX = TWO/EX
              NCALC = 2
              DO 90 I = 3,NB
                  IF (TWOBYX.LT.ONE) THEN
                      IF (ABS(BY(I-1))*TWOBYX.GE.XINF/AYE) GO TO 100

                  ELSE
                      IF (ABS(BY(I-1)).GE.XINF/AYE/TWOBYX) GO TO 100
                  END IF

                  BY(I) = TWOBYX*AYE*BY(I-1) - BY(I-2)
                  AYE = AYE + ONE
                  NCALC = NCALC + 1
   90         CONTINUE
          END IF

  100     DO 110 I = NCALC + 1,NB
              BY(I) = ZERO
  110     CONTINUE

      ELSE
          BY(1) = ZERO
          NCALC = MIN(NB,0) - 1
      END IF

      RETURN
C---------- Last line of RYBESL ----------
      END
C     END OF SUBROUTINE *** RYBESL ***
C
C     BEGIN OF FUNCTION *** DLGAMA ***
CS    REAL FUNCTION ALGAMA(X)
      DOUBLE PRECISION FUNCTION DLGAMA(X)
C
C     ##################################################################
C
C     *** T. WIEDER, 05.03.1998 ***
C
C     SOURCE FOR *** FUNCTION DLGAMA ***:
C
C     SUBPROGRAM LIBRARY SPECFUN
C
C     AVAILABILITIY: PUBLIC DOMAIN
C
C     DEVELOPER: W.J. CODY AND L. STOLTZ, APPLED MATHEMATICS DIVISION,
C     ARGONNE NATIONAL LABORATORY, ARGONNE, IL 60439.
C
C
C     DISTRIBUTOR: NETLIB
C
C     ##################################################################
C
C----------------------------------------------------------------------
C
C This routine calculates the LOG(GAMMA) function for a positive real
C   argument X.  Computation is based on an algorithm outlined in
C   references 1 and 2.  The program uses rational functions that
C   theoretically approximate LOG(GAMMA) to at least 18 significant
C   decimal digits.  The approximation for X > 12 is from reference
C   3, while approximations for X < 12.0 are similar to those in
C   reference 1, but are unpublished.  The accuracy achieved depends
C   on the arithmetic system, the compiler, the intrinsic functions,
C   and proper selection of the machine-dependent constants.
C
C
C*********************************************************************
C*********************************************************************
C
C Explanation of machine-dependent constants
C
C beta   - radix for the floating-point representation
C maxexp - the smallest positive power of beta that overflows
C XBIG   - largest argument for which LN(GAMMA(X)) is representable
C          in the machine, i.e., the solution to the equation
C                  LN(GAMMA(XBIG)) = beta**maxexp
C XINF   - largest machine representable floating-point number;
C          approximately beta**maxexp.
C EPS    - The smallest positive floating-point number such that
C          1.0+EPS .GT. 1.0
C FRTBIG - Rough estimate of the fourth root of XBIG
C
C
C     Approximate values for some important machines are:
C
C                           beta      maxexp         XBIG
C
C CRAY-1        (S.P.)        2        8191       9.62E+2461
C Cyber 180/855
C   under NOS   (S.P.)        2        1070       1.72E+319
C IEEE (IBM/XT,
C   SUN, etc.)  (S.P.)        2         128       4.08E+36
C IEEE (IBM/XT,
C   SUN, etc.)  (D.P.)        2        1024       2.55D+305
C IBM 3033      (D.P.)       16          63       4.29D+73
C VAX D-Format  (D.P.)        2         127       2.05D+36
C VAX G-Format  (D.P.)        2        1023       1.28D+305
C
C
C                           XINF        EPS        FRTBIG
C
C CRAY-1        (S.P.)   5.45E+2465   7.11E-15    3.13E+615
C Cyber 180/855
C   under NOS   (S.P.)   1.26E+322    3.55E-15    6.44E+79
C IEEE (IBM/XT,
C   SUN, etc.)  (S.P.)   3.40E+38     1.19E-7     1.42E+9
C IEEE (IBM/XT,
C   SUN, etc.)  (D.P.)   1.79D+308    2.22D-16    2.25D+76
C IBM 3033      (D.P.)   7.23D+75     2.22D-16    2.56D+18
C VAX D-Format  (D.P.)   1.70D+38     1.39D-17    1.20D+9
C VAX G-Format  (D.P.)   8.98D+307    1.11D-16    1.89D+76
C
C**************************************************************
C**************************************************************
C
C Error returns
C
C  The program returns the value XINF for X .LE. 0.0 or when
C     overflow would occur.  The computation is believed to
C     be free of underflow and overflow.
C
C
C Intrinsic functions required are:
C
C      LOG
C
C
C References:
C
C  1) W. J. Cody and K. E. Hillstrom, 'Chebyshev Approximations for
C     the Natural Logarithm of the Gamma Function,' Math. Comp. 21,
C     1967, pp. 198-203.
C
C  2) K. E. Hillstrom, ANL/AMD Program ANLC366S, DGAMMA/DLGAMA, May,
C     1969.
C
C  3) Hart, Et. Al., Computer Approximations, Wiley and sons, New
C     York, 1968.
C
C
C  Authors: W. J. Cody and L. Stoltz
C           Argonne National Laboratory
C
C  Latest modification: June 16, 1988
C
C----------------------------------------------------------------------
CS    REAL
C     .. Scalar Arguments ..
      DOUBLE PRECISION X
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION CORR,D1,D2,D4,EPS,FOUR,FRTBIG,HALF,ONE,PNT68,RES,
     +                 SQRTPI,THRHAL,TWELVE,TWO,XBIG,XDEN,XINF,XM1,XM2,
     +                 XM4,XNUM,Y,YSQ,ZERO
      INTEGER I
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION C(7),P1(8),P2(8),P4(8),Q1(8),Q2(8),Q4(8)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC LOG
C     ..
C     .. Data statements ..
C----------------------------------------------------------------------
C  Mathematical constants
C----------------------------------------------------------------------
CS    DATA ONE,HALF,TWELVE,ZERO/1.0E0,0.5E0,12.0E0,0.0E0/,
CS   1     FOUR,THRHAL,TWO,PNT68/4.0E0,1.5E0,2.0E0,0.6796875E0/,
CS   2     SQRTPI/0.9189385332046727417803297E0/
C----------------------------------------------------------------------
C  Machine dependent parameters
C----------------------------------------------------------------------
CS    DATA XBIG,XINF,EPS,FRTBIG/4.08E36,3.401E38,1.19E-7,1.42E9/
C----------------------------------------------------------------------
C  Numerator and denominator coefficients for rational minimax
C     approximation over (0.5,1.5).
C----------------------------------------------------------------------
CS    DATA D1/-5.772156649015328605195174E-1/
CS    DATA P1/4.945235359296727046734888E0,2.018112620856775083915565E2,
CS   1        2.290838373831346393026739E3,1.131967205903380828685045E4,
CS   2        2.855724635671635335736389E4,3.848496228443793359990269E4,
CS   3        2.637748787624195437963534E4,7.225813979700288197698961E3/
CS    DATA Q1/6.748212550303777196073036E1,1.113332393857199323513008E3,
CS   1        7.738757056935398733233834E3,2.763987074403340708898585E4,
CS   2        5.499310206226157329794414E4,6.161122180066002127833352E4,
CS   3        3.635127591501940507276287E4,8.785536302431013170870835E3/
C----------------------------------------------------------------------
C  Numerator and denominator coefficients for rational minimax
C     Approximation over (1.5,4.0).
C----------------------------------------------------------------------
CS    DATA D2/4.227843350984671393993777E-1/
CS    DATA P2/4.974607845568932035012064E0,5.424138599891070494101986E2,
CS   1        1.550693864978364947665077E4,1.847932904445632425417223E5,
CS   2        1.088204769468828767498470E6,3.338152967987029735917223E6,
CS   3        5.106661678927352456275255E6,3.074109054850539556250927E6/
CS    DATA Q2/1.830328399370592604055942E2,7.765049321445005871323047E3,
CS   1        1.331903827966074194402448E5,1.136705821321969608938755E6,
CS   2        5.267964117437946917577538E6,1.346701454311101692290052E7,
CS   3        1.782736530353274213975932E7,9.533095591844353613395747E6/
C----------------------------------------------------------------------
C  Numerator and denominator coefficients for rational minimax
C     Approximation over (4.0,12.0).
C----------------------------------------------------------------------
CS    DATA D4/1.791759469228055000094023E0/
CS    DATA P4/1.474502166059939948905062E4,2.426813369486704502836312E6,
CS   1        1.214755574045093227939592E8,2.663432449630976949898078E9,
CS   2      2.940378956634553899906876E10,1.702665737765398868392998E11,
CS   3      4.926125793377430887588120E11,5.606251856223951465078242E11/
CS    DATA Q4/2.690530175870899333379843E3,6.393885654300092398984238E5,
CS   2        4.135599930241388052042842E7,1.120872109616147941376570E9,
CS   3      1.488613728678813811542398E10,1.016803586272438228077304E11,
CS   4      3.417476345507377132798597E11,4.463158187419713286462081E11/
C----------------------------------------------------------------------
C  Coefficients for minimax approximation over (12, INF).
C----------------------------------------------------------------------
CS    DATA C/-1.910444077728E-03,8.4171387781295E-04,
CS   1     -5.952379913043012E-04,7.93650793500350248E-04,
CS   2     -2.777777777777681622553E-03,8.333333333333333331554247E-02,
CS   3      5.7083835261E-03/
      DATA ONE,HALF,TWELVE,ZERO/1.0D0,0.5D0,12.0D0,0.0D0/,FOUR,THRHAL,
     +     TWO,PNT68/4.0D0,1.5D0,2.0D0,0.6796875D0/,
     +     SQRTPI/0.9189385332046727417803297D0/
      DATA XBIG,XINF,EPS,FRTBIG/2.55D305,1.79D308,2.22D-16,2.25D76/
      DATA D1/-5.772156649015328605195174D-1/
      DATA P1/4.945235359296727046734888D0,2.018112620856775083915565D2,
     +     2.290838373831346393026739D3,1.131967205903380828685045D4,
     +     2.855724635671635335736389D4,3.848496228443793359990269D4,
     +     2.637748787624195437963534D4,7.225813979700288197698961D3/
      DATA Q1/6.748212550303777196073036D1,1.113332393857199323513008D3,
     +     7.738757056935398733233834D3,2.763987074403340708898585D4,
     +     5.499310206226157329794414D4,6.161122180066002127833352D4,
     +     3.635127591501940507276287D4,8.785536302431013170870835D3/
      DATA D2/4.227843350984671393993777D-1/
      DATA P2/4.974607845568932035012064D0,5.424138599891070494101986D2,
     +     1.550693864978364947665077D4,1.847932904445632425417223D5,
     +     1.088204769468828767498470D6,3.338152967987029735917223D6,
     +     5.106661678927352456275255D6,3.074109054850539556250927D6/
      DATA Q2/1.830328399370592604055942D2,7.765049321445005871323047D3,
     +     1.331903827966074194402448D5,1.136705821321969608938755D6,
     +     5.267964117437946917577538D6,1.346701454311101692290052D7,
     +     1.782736530353274213975932D7,9.533095591844353613395747D6/
      DATA D4/1.791759469228055000094023D0/
      DATA P4/1.474502166059939948905062D4,2.426813369486704502836312D6,
     +     1.214755574045093227939592D8,2.663432449630976949898078D9,
     +     2.940378956634553899906876D10,1.702665737765398868392998D11,
     +     4.926125793377430887588120D11,5.606251856223951465078242D11/
      DATA Q4/2.690530175870899333379843D3,6.393885654300092398984238D5,
     +     4.135599930241388052042842D7,1.120872109616147941376570D9,
     +     1.488613728678813811542398D10,1.016803586272438228077304D11,
     +     3.417476345507377132798597D11,4.463158187419713286462081D11/
      DATA C/-1.910444077728D-03,8.4171387781295D-04,
     +     -5.952379913043012D-04,7.93650793500350248D-04,
     +     -2.777777777777681622553D-03,8.333333333333333331554247D-02,
     +     5.7083835261D-03/
C     ..
C----------------------------------------------------------------------
      Y = X
      IF ((Y.GT.ZERO) .AND. (Y.LE.XBIG)) THEN
          IF (Y.LE.EPS) THEN
              RES = -LOG(Y)

          ELSE IF (Y.LE.THRHAL) THEN
C----------------------------------------------------------------------
C  EPS .LT. X .LE. 1.5
C----------------------------------------------------------------------
              IF (Y.LT.PNT68) THEN
                  CORR = -LOG(Y)
                  XM1 = Y

              ELSE
                  CORR = ZERO
                  XM1 = (Y-HALF) - HALF
              END IF

              IF ((Y.LE.HALF) .OR. (Y.GE.PNT68)) THEN
                  XDEN = ONE
                  XNUM = ZERO
                  DO 10 I = 1,8
                      XNUM = XNUM*XM1 + P1(I)
                      XDEN = XDEN*XM1 + Q1(I)
   10             CONTINUE
                  RES = CORR + (XM1* (D1+XM1* (XNUM/XDEN)))

              ELSE
                  XM2 = (Y-HALF) - HALF
                  XDEN = ONE
                  XNUM = ZERO
                  DO 20 I = 1,8
                      XNUM = XNUM*XM2 + P2(I)
                      XDEN = XDEN*XM2 + Q2(I)
   20             CONTINUE
                  RES = CORR + XM2* (D2+XM2* (XNUM/XDEN))
              END IF

          ELSE IF (Y.LE.FOUR) THEN
C----------------------------------------------------------------------
C  1.5 .LT. X .LE. 4.0
C----------------------------------------------------------------------
              XM2 = Y - TWO
              XDEN = ONE
              XNUM = ZERO
              DO 30 I = 1,8
                  XNUM = XNUM*XM2 + P2(I)
                  XDEN = XDEN*XM2 + Q2(I)
   30         CONTINUE
              RES = XM2* (D2+XM2* (XNUM/XDEN))

          ELSE IF (Y.LE.TWELVE) THEN
C----------------------------------------------------------------------
C  4.0 .LT. X .LE. 12.0
C----------------------------------------------------------------------
              XM4 = Y - FOUR
              XDEN = -ONE
              XNUM = ZERO
              DO 40 I = 1,8
                  XNUM = XNUM*XM4 + P4(I)
                  XDEN = XDEN*XM4 + Q4(I)
   40         CONTINUE
              RES = D4 + XM4* (XNUM/XDEN)

          ELSE
C----------------------------------------------------------------------
C  Evaluate for argument .GE. 12.0,
C----------------------------------------------------------------------
              RES = ZERO
              IF (Y.LE.FRTBIG) THEN
                  RES = C(7)
                  YSQ = Y*Y
                  DO 50 I = 1,6
                      RES = RES/YSQ + C(I)
   50             CONTINUE
              END IF

              RES = RES/Y
              CORR = LOG(Y)
              RES = RES + SQRTPI - HALF*CORR
              RES = RES + Y* (CORR-ONE)
          END IF

      ELSE
C----------------------------------------------------------------------
C  Return for bad arguments
C----------------------------------------------------------------------
          RES = XINF
      END IF
C----------------------------------------------------------------------
C  Final adjustments and return
C----------------------------------------------------------------------
CS    ALGAMA = RES
      DLGAMA = RES
      RETURN
C ---------- Last line of DLGAMA ----------
      END
C     END OF FUNCTION *** DLGAMA ***
C
C     BEGIN OF SUBROUTINE *** INTHP ***
      SUBROUTINE INTHP(A,B,D,F,M,P,EPS,INF,QUADR)
C
C     ##################################################################
C
C     *** T. WIEDER, 05.03.1998 ***
C
C     SOURCE FOR SUBROUTINE *** INTHP ***:
C
C     SUBPROGRAM LIBRARY TOMS
C
C     AVAILABILITIY: PUBLIC DOMAIN
C
C     DEVELOPER: K. SIKORSKY, F. STENGER, AND J. SCHWING.
C
C     DISTRIBUTOR: TRANSACTIONS ON MATHEMATICAL SOFTWARE,
C     VOL. 10 (1984), PP. 152 - 160.
C
C     ##################################################################
C
C     ALGORITHM 614 COLLECTED ALGORITHMS FROM ACM.
C     ALGORITHM APPEARED IN ACM-TRANS. MATH. SOFTWARE, VOL.10, NO. 2,
C     JUN., 1984, P. 152-160.
C
C        THIS SUBROUTINE COMPUTES INTEGRAL OF FUNCTIONS WHICH
C     MAY HAVE SINGULARITIES AT ONE OR BOTH END-POINTS OF AN
C     INTERVAL (A,B), SEE [1, 2]. IT CONTAINS FOUR DIFFERENT
C     QUADRATURE ROUTINES: ONE OVER A FINITE INTERVAL (A,B),
C     TWO OVER (A,+INFINITY), AND ONE OVER (-INFINITY,+INFINITY).
C     OF THE TWO FORMULAS OVER (A,+INFINITY), THE FIRST (INF=2
C     BELOW) IS MORE SUITED TO NON-OSCILLATORY INTEGRANDS, WHILE
C     THE SECOND (INF=3) IS MORE SUITED TO OSCILLATORY INTEGRANDS.
C        THE USER SUPPLIES THE INTEGRAND FUNCTION, HE SPECIFIES THE
C     INTERVAL, AS WELL AS THE RELATIVE ERROR TO WHICH THE INTE-
C     GRAL IS TO BE EVALUATED.
C        THE FORMULAS ARE OPTIMAL IN CERTAIN HARDY SPACES H(P,DD),
C     SEE [1, 2]. HERE DD IS AN OPEN DOMAIN IN THE COMPLEX PLANE,
C     A AND B BELONG TO THE BOUNDARY OF DD AND H(P,DD), P.GT.1, IS
C     THE SET OF ALL ANALYTIC FUNCTONS IN DD WHOSE P-TH NORM DEFI-
C     NED AS IN [2] IS FINITE.
C        IF THE USER IS UNABLE TO SPECIFY THE PARAMETERS P AND D
C     OF THE SPACE H(P,DD) TO WHICH HIS INTEGRAND BELONGS, THE
C     ALGORITHM TERMINATES ACCORDING TO A HEURISTIC CRITERION, SEE
C     [2] AND COMMENTS TO EPS.
C        IF THE USER CAN SPECIFY THE PARAMETERS P AND D OF THE
C     SPACE H(P,DD) TO WHICH HIS INTEGRAND BELONGS, THE ALGORITHM
C     TERMINATES WITH AN ANSWER HAVING A GUARANTEED ACCURACY ( DE-
C     TEMINISTIC CRITERION, SEE [1, 2] AND COMMENTS TO EPS).
C
C
C     INPUT PARAMETERS
C
C
C     A = LOWER LIMIT OF INTEGRATION (SEE COMMENTS TO INF).
C
C     B = UPPER LIMIT OF INTEGRATION (SEE COMMENTS TO INF).
C
C     D = A PARAMETER OF THE CLASS H(P,DD) (SEE COMMENTS TO
C         INF).
C
C         USER SETS D:
C
C         HEURISTIC TERMINATION
C       = ANY REAL NUMBER
C
C         DETERMINISTIC TERMINATION
C       = A NUMBER IN THE RANGE 0.LT.D.LE.PI/2.
C
C     F = A NAME OF AN EXTERNAL INTEGRAND FUNCTION TO BE
C         SUPPLIED BY THE USER. F(X) COMPUTES THE VALUE OF
C         A FUNCTION F AT A POINT X. THE STATEMENT
C         ...EXTERNAL F... MUST APPEAR IN THE MAIN PROGRAM.
C
C     M = MAXIMAL NUMBER OF FUNCTION EVALUATIONS ALLOWED IN
C         THE COMPUTATIONS, M.GE.3.( ALTERED ON EXIT ).
C
C     P = 0, 1, .GT.1  A PARAMETER OF THE CLASS H(P,DD).
C
C         USER SETS P:
C       = 0 - HEURISTIC TERMINATION.
C       = 1 - DETERMINISTIC TERMINATION WITH THE INFINITY
C             NORM.
C      .GT.1 -DETERMINISTIC TERMINATION WITH THE P-TH NORM.
C
C   EPS = A REAL NUMBER - THE RELATIVE ERROR BOUND - SEE
C         REMARKS BELOW. ( ALTERED ON EXIT ).
C
C
C   INF = 1, 2, 3, 4 - INFORMATION PARAMETER. ( ALTERED ON EXIT ).
C
C       = 1 SIGNIFIES AN INFINITE INTERVAL (A,B)=REAL LINE,
C           A AND B ANY NUMBERS.
C           DETERMINISTIC TERMINATION -
C           DD=STRIP(Z:ABS(IM(Z)).LT.D).
C
C       = 2 SIGNIFIES A SEMI-INFINITE INTERVAL (A, +INFINITY)
C           USER SUPPLIES A, B ANY NUMBER.
C           QUADRATURE SUITED TO NON-OSCILLATORY INTEGRANDS.
C           DETERMINISTIC TERMINATION -
C           DD=SECTOR(Z:ABS(ARG(Z-A)).LT.D).
C
C       = 3 SIGNIFIES A SEMI INFINITE INTERVAL (A,+INFINITY)
C           USER SUPPLIES A, B ANY NUMBER.
C           QUADRATURE SUITED TO OSCILLATORY INTEGRANDS.
C           DETERMINISTIC TERMINATION -
C           DD=REGION(Z:ABS(ARG(SINH(Z-A))).LT.D).
C
C       = 4 SIGNIFIES A FINITE INTERVAL (A,B).
C           USER SUPPLIES A AND B.
C           DETERMINISTIC TERMINATION -
C           DD=LENS REGION(Z:ABS(ARG((Z-A)/(B-Z))).LT.D).
C
C
C     OUTPUT PARAMETERS
C
C
C     M = THE NUMBER OF FUNCTION EVALUATIONS USED IN THE
C         QUADRATURE.
C
C   EPS = THE RELATIVE ERROR BOUND (SEE REMARKS BELOW).
C
C         DETERMINISTIC TERMINATION
C
C       = THE RELATIVE ERROR REXA BOUND, I.E.,
C                 REXA(F,M(OUTPUT)) .LE. EPS.
C
C         HEURISTIC TERMINATION
C
C       = MAX(EPS(INPUT),MACHEP).
C
C   INF = 0, 1 - DETERMINISTIC TERMINATION
C
C       = 0 COMPUTED QUADRATURE QCOM(F,M(EPS)), SEE REMARKS
C           BELOW.
C
C       = 1 COMPUTED QUADRATURE QCOM(F,M1), SEE REMARKS
C           BELOW.
C
C   INF = 2, 3, 4 - HEURISTIC TERMINATION.
C
C       = 2 INTEGRATION COMPLETED WITH EPS=MAX(EPS(INPUT),
C           MACHEP). WE CAN EXPECT THE RELATIVE ERROR
C           REXA TO BE OF THE ORDER OF EPS (FOR SOME P.GE.1).
C
C       = 3 INTEGRATION NOT COMPLETED. ATTEMPT TO EXCEED THE
C           MAXIMAL ALLOWED NUMBER OF FUNCTION EVALUATIONS M.
C           TRUNCATION CONDITIONS (SEE [2]) SATISFIED. QUADR
C           SET TO BE EQUAL TO THE LAST TRAPEZOIDAL APPRO-
C           XIMATION. IT IS LIKELY THAT QUADR APPROXIMATES THE
C           INTEGRAL IF M IS LARGE.
C
C       = 4 INTEGRATION NOT COMPLETED. ATTEMPT TO EXCEED THE
C           MAXIMAL ALLOWED NUMBER OF FUNCTION EVALUATIONS M.
C           TRUNCATION CONDITIONS (SEE [2]) NOT SATISFIED.
C           QUADR SET TO BE EQUAL TO THE COMPUTED TRAPEZOIDAL
C           APPROXIMATION. IT IS UNLIKELY THAT QUADR APPROXIMATES
C           THE INTEGRAL.
C
C   INF = 10, 11, 12, 13 - INCORRECT INPUT
C
C       = 10  M.LT.3.
C
C       = 11  P DOES NOT SATISFY P=0, P=1 OR P.GT.1 OR IN THE
C             CASE OF DETERMINISTIC TERMINATION D DOES NOT
C             SATISFY 0.LT.D.LE.PI/2.
C
C       = 12  A.GE.B IN CASE OF A FINITE INTERVAL.
C
C       = 13  INF NOT EQUAL TO 1, 2, 3, OR 4.
C
C
C   QUADR = THE COMPUTED VALUE OF QUADRATURE.
C
C
C     REMARKS:
C
C         LET  QEXA(F,M)  ( QCOM(F,M) ) BE THE EXACT (COMPUTED)
C         VALUE OF THE QUADRATURE WITH M FUNCTION EVALUATIONS,
C         AND LET  REXA(F,M) ( RCOM(F,M) ) BE THE RELATIVE ERROR
C         OF QEXA (QCOM) ,I.E.,
C            REXA(F,M)=ABS(INTEGRAL(F)-QEXA(F,M))/NORM(F),
C            RCOM(F,M)=ABS(INTEGRAL(F)-QCOM(F,M))/NORM(F),
C         WITH THE NOTATION 0/0=0.
C             DUE TO THE ROUNDOFF ONE CANNOT EXPECT THE ERROR
C         RCOM TO BE LESS THAN THE RELATIVE MACHINE PRECISION
C         MACHEP. THEREFORE THE INPUT VALUE OF EPS IS CHANGED
C         ACCORDING TO THE FORMULA
C                   EPS=MAX(EPS,MACHEP).
C
C         DETERMINISTIC TERMINATON CASE
C
C             THE NUMBER OF FUNCTON EVALUATIONS M(EPS) IS COMPUTED
C         SO THAT THE ERROR REXA IS NO GREATER THAN EPS,I.E.,
C
C         (*)     REXA(F,M(EPS)) .LE. EPS .
C
C         IF M(EPS).LE.M THEN THE QUADRATURE QCOM(F,M(EPS)) IS COM-
C         PUTED. OTHERWISE, WHICH MEANS THAT EPS IS TOO SMALL WITH
C         RESPECT TO M, THE QUADRATURE QCOM(F,M1) IS COMPUTED, WHERE
C         M1=2*INT((M-1)/2)+1. IN THIS CASE EPS IS CHANGED TO THE
C         SMALLEST NUMBER FOR WHICH THE ESTIMATE (*) HOLDS WITH
C         M(EPS)=M1 FUNCTION EVALUATIONS.
C
C         HEURISTIC TERMINATION CASE
C
C             WE CAN EXPECT THE RELATIVE ERROR REXA TO BE OF THE
C         ORDER OF EPS, SEE [2]. IF EPS IS TOO SMALL WITH RESPECT
C         TO M THEN THE QUADRATURE QCOM(F,M) IS COMPUTED.
C
C         ROUNDOFF ERRORS
C
C             IN BOTH DETERMINISTIC AND HEURISTIC CASES THE ROUND-
C         OFF ERROR
C                    ROFF=ABS(QEXA(F,M)-QCOM(F,M))
C         CAN BE ESTIMATED BY
C
C         (**)       ROFF .LE. 3*C1*R*MACHEP,
C
C         WHERE  R=QCOM(ABS(F),M)+(1+2*C2)/3*SUM(W(I),I=1,2,...M)
C         AND C1 IS OF THE ORDER OF UNITY, C1=1/(1-3*MACHEP), W(I)
C         ARE THE WEIGHTS OF THE QUADRATURE, SEE [2], AND C2 IS
C         A CONSTANT ESTIMATING THE ACCURACY OF COMPUTING FUNCTION
C         VALUES, I.E.,
C               ABS(EXACT(F(X))-COMPUTED(F(X))).LE.C2*MACHEP.
C         IF THE INTEGRAND VALUES ARE COMPUTED INACCURATELY, I.E.,
C         C2 IS LARGE, THEN THE ESTIMATE (**) IS LARGE AND ONE CAN
C         EXPECT THE ACTUAL ERROR ROFF TO BE LARGE. NUMERICAL TESTS
C         INDICATE THAT THIS HAPPENS ESPECIALLY WHEN THE INTEGRAND
C         IS EVALUATED INACCURATELY NEAR A SINGULARITY. THE WAYS OF
C         CIRCUMVENTING SUCH PITFALLS ARE EXPLAINED IN [2].
C
C     REFERENCES:
C
C     [1] SIKORSKI,K., OPTIMAL QUADRATURE ALGORITHMS IN HP
C            SPACES, NUM. MATH., 39, 405-410 (1982).
C     [2] SIKORSKI,K., STENGER,F., OPTIMAL QUADRATURES IN
C            HP SPACES, ACM TOMS.
C
C
C      IMPLICIT NONE
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION A,B,D,EPS,P,QUADR
      INTEGER INF,M
C     ..
C     .. Function Arguments ..
      DOUBLE PRECISION F
      EXTERNAL F
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ALFA,BA,C,C0,COR,E1,EPS3,EXPH,EXPH0,H,H0,H1,PI,S,
     +                 S1,SQ2,SR,SUM,SUM1,SUM2,T,U,V,V0,V1,V2,W,W1,W2,
     +                 W3,W4
      INTEGER I,I1,K,L,L1,M1,M2,N,N1
      LOGICAL INF1,INF2
C     ..
C     .. External Functions ..
      DOUBLE PRECISION D1MACH
      EXTERNAL D1MACH
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,DBLE,DLOG,DSQRT,DTAN,EXP,FLOAT,INT,SQRT
C     ..
      PI = 4.0D0*DTAN(1.0D0)
C
C     CHECK THE INPUT DATA
C
      IF (INF.NE.1 .AND. INF.NE.2 .AND. INF.NE.3 .AND.
     +    INF.NE.4) GO TO 300
      IF (M.LT.3) GO TO 270
      IF (P.LT.1.0D0 .AND. P.NE.0.0D0) GO TO 280
      IF (P.GE.1.0D0 .AND. (D.LE.0.D0.OR.D.GT.PI/2.D0)) GO TO 280
      IF (INF.EQ.4 .AND. A.GE.B) GO TO 290
C
      SQ2 = DSQRT(2.0D0)
      I1 = INF - 2
      BA = B - A
      N1 = 0
C
C     COMPUTE THE RELATIVE MACHINE PRECISION AND CHECK
C     THE VALUE OF EPS.  CAUTION...THIS LOOP MAY NOT WORK ON A
C     MACHINE THAT HAS AN ACCURATED ARITHMETIC PROCESS COMPARED
C     TO THE STORAGE PRECISION.  THE VALUE OF U MAY NEED TO BE
C     SIMPLY DEFINED AS THE RELATIVE ACCURACY OF STORAGE PRECISION.
C
      U = 1.0D0
   10 U = U/10.0D0
      T = 1.0D0 + U
      IF (1.0D0.NE.T) GO TO 10
      U = U*10.0D0
      IF (EPS.LT.U) EPS = U
C
      IF (P.EQ.0.D0) GO TO 40
C
C     SET UP DATA FOR THE DETERMINISTIC TERMINATION
C
      IF (P.EQ.1.0D0) ALFA = 1.0D0
      IF (P.GT.1.0D0) ALFA = (P-1.0D0)/P
      C = 2.0D0*PI/ (1.0D0-1.D0/EXP(PI*DSQRT(ALFA))) + 4.0D0**ALFA/ALFA
      W = DLOG(C/EPS)
      W1 = 1.0D0/ (PI*PI*ALFA)*W*W
      N = INT(W1)
      IF (W1.GT.FLOAT(N)) N = N + 1
      IF (W1.EQ.0.D0) N = 1
      N1 = 2*N + 1
      SR = DSQRT(ALFA*DBLE(N))
      IF (N1.LE.M) GO TO 20
C
C     EPS TOO SMALL WITH RESPECT TO M. COMPUTE THE NEW EPS
C     GUARANTEED BY THE VALUE OF M.
C
      N1 = 1
      N = INT(FLOAT((M-1)/2))
      SR = DSQRT(ALFA*DBLE(N))
      M = 2*N + 1
      EPS = C/EXP(PI*SR)
      GO TO 30
C
   20 M = N1
      N1 = 0
   30 H = 2.0D0*D/SR
      SUM2 = 0.0D0
      L1 = N
      K = N
      INF1 = .FALSE.
      INF2 = .FALSE.
      H0 = H
      GO TO 50
C
C     SET UP DATA FOR THE HEURISTIC TERMINATION
C
   40 H = 1.0D0
      H0 = 1.0D0
      EPS3 = EPS/3.0D0
      SR = SQRT(EPS)
      V1 = EPS*10.0D0
      V2 = V1
      M1 = M - 1
      N = INT(FLOAT(M1/2))
      M2 = N
      L1 = 0
      INF1 = .TRUE.
      INF2 = .FALSE.
C
C     INITIALIZE THE QUADRATURE
C
   50 I = 0
      IF (INF.EQ.1) SUM = F(0.0D0)
      IF (INF.EQ.2) SUM = F(A+1.0D0)
      IF (INF.EQ.3) SUM = F(A+DLOG(1.0D0+SQ2))/SQ2
      IF (INF.EQ.4) SUM = F((A+B)/2.0D0)/4.0D0*BA
C
C     COMPUTE WEIGHTS, NODES AND FUNCTION VALUES
C
   60 EXPH = EXP(H)
      EXPH0 = EXP(H0)
      H1 = H0
      E1 = EXPH0
      U = 0.0D0
      COR = 0.0D0
C
C     *** T. WIEDER, 24.02.1998 ***
C70    IF (I1) 80,90,100
   70 IF (I1.LT.0) THEN
          GO TO 80

      ELSE IF (I1.EQ.0) THEN
          GO TO 90

      ELSE IF (I1.GT.0) THEN
          GO TO 100

      END IF
C
   80 V = F(H1)
      H1 = H1 + H
      GO TO 150
C
   90 V = E1*F(A+E1)
      E1 = E1*EXPH
      GO TO 150
C
  100 IF (INF.EQ.4) GO TO 140
      W1 = DSQRT(E1+1.0D0/E1)
      W2 = DSQRT(E1)
      IF (E1.LT.0.1D0) GO TO 110
      S = DLOG(E1+W1*W2)
      GO TO 130

  110 W3 = E1
      W4 = E1*E1
      C0 = 1.0D0
      S = E1
      S1 = E1
      T = 0.0D0
  120 C0 = -C0* (0.5D0+T)* (2.0D0*T+1.D0)/ (2.0D0*T+3.0D0)/ (T+1.0D0)
      T = T + 1.0D0
      W3 = W3*W4
      S = S + C0*W3
      IF (S.EQ.S1) GO TO 130
      S1 = S
      GO TO 120

  130 V = W2/W1*F(A+S)
C     *** 05.02.1998 ***
      IF (ABS(E1).GT.D1MACH(2)/10.0D0) THEN
          E1 = E1/10.0D0
          E1 = E1*EXPH

      ELSE
          E1 = E1*EXPH
      END IF
C     *** 05.02.1998 ***
      GO TO 150
C
  140 W1 = E1 + 1.0D0
      V = E1/W1/W1*F((A+B*E1)/W1)*BA
      E1 = E1*EXPH
C
C     SUMMATION ALGORITHM
C
  150 I = I + 1
      SUM1 = U + V
      IF (ABS(U).LT.ABS(V)) GO TO 160
      COR = V - (SUM1-U) + COR
      GO TO 170

  160 COR = U - (SUM1-V) + COR
  170 U = SUM1
      IF (I.LT.L1) GO TO 70
C
C     SWITCH TO CHECK TRUNCATION CONDITION ( HEURISTIC
C     TERMINATION)
C
      IF (INF1) GO TO 190
C
C     SWITCH TO COMPUTE THE MIDORDINATE APPROXIMATION
C     ( HEURISTIC TERMINATION ) OR TO STOP ( DETERMINIS-
C     TIC TERMINATION)
C
      IF (INF2) GO TO 210
C
C     SET UP PARAMETERS TO CONTINUE SUMMATION
C
      L1 = K
  180 INF2 = .TRUE.
      I = 0
      EXPH = 1.D0/EXPH
      H0 = -H0
      E1 = 1.0D0/EXPH0
      H1 = H0
      H = -H
      GO TO 70
C
C     TRUNCATION CONDITION
C
  190 V0 = V1
      V1 = V2
      V2 = ABS(V)
      IF (V0+V1+V2.LE.EPS3) GO TO 200
      IF (I.LT.M2) GO TO 70
      N1 = 5
  200 IF (INF2) K = I
      IF (.NOT.INF2) L = I
      V1 = 10.D0*EPS
      V2 = V1
      M2 = M1 - L
      IF (.NOT.INF2) GO TO 180
C
C     N1=5 - TRUNCATION CONDITION NOT SATISFIED
C
      IF (N1.EQ.5) GO TO 260
C
C     TRUNCATION CONDITION SATISFIED, SUM2=TRAPEZOIDAL
C     APPROXIMATION
C
      SUM2 = SUM1 + COR + SUM
      M2 = 2* (K+L)
C
C     CHECK THE NUMBER OF FUNCTION EVALUATIONS
C
      IF (M2.GT.M1) GO TO 240
C
C     INITIALIZE ITERATION
C
      INF1 = .FALSE.
      INF2 = .FALSE.
      L1 = L
      I = 0
      H = -H
      H0 = H/2.D0
      GO TO 60
C
C     P.GE.1 = DETERMINISTIC TERMINATION
C
  210 IF (P.GE.1.D0) GO TO 220
C
C     COMPUTE THE MIDORDINATE APPROXIMATION SUM1
C
      H = -H
      SUM1 = (SUM1+COR)*H
      W1 = (SUM1+SUM2)/2.0D0
C
C     TERMINATION CONDITION
C
      IF (ABS(SUM1-SUM2).LE.SR) GO TO 230
C
C     SET UP DATA FOR THE NEXT ITERATION
C
      M2 = 2*M2
      IF (M2.GT.M1) GO TO 250
      I = 0
      K = 2*K
      L = 2*L
      L1 = L
      H = H/2.0D0
      H0 = H/2.0D0
      SUM2 = W1
      INF2 = .FALSE.
      GO TO 60
C
C     FINAL RESULTS
C
  220 QUADR = -H* (SUM1+COR+SUM)
      INF = N1
      RETURN
C
  230 QUADR = W1
      INF = 2
      M = M2 + 1
      RETURN
C
  240 QUADR = SUM2
      INF = 3
      M = K + L + 1
      RETURN
C
  250 QUADR = W1
      INF = 3
      M = M2/2 + 1
      RETURN
C
  260 QUADR = U + COR + SUM
      INF = 4
      M = K + L + 1
      RETURN
C
  270 INF = 10
      RETURN
C
  280 INF = 11
      RETURN
C
  290 INF = 12
      RETURN
C
  300 INF = 13
      RETURN
C
      END
C     END OF SUBROUTINE *** INTHP ***
C
C     ##################################################################
C     ##################################################################
C
C     BEGIN OF SUBROUTINE *** OSCINT ***
      SUBROUTINE OSCINT(AZERO,PERIOD,RFIRST,EPS,NQUAD,NDIM1,NDIM2,GAUSS,
     +                  HFUN,GPER,WORK,SAVPER,WEIGHT,ABSCIS,QLIST,
     +                  RESULT,ISTATE)
C
C     ##################################################################
C
C     *** T. WIEDER, 05.03.1998 ***
C
C     SOURCE FOR SUBROUTINE *** OSCINT ***:
C
C     SUBPROGRAM LIBRARY TOMS
C
C     AVAILABILITIY: PUBLIC DOMAIN
C
C     DEVELOPER: JAMES LYNESS AND GWENDOLEN HINES.
C
C     DISTRIBUTOR: TRANSACTIONS ON MATHEMATICAL SOFTWARE,
C     VOL. 12 (1986), PP. 24 - 25.
C
C     ##################################################################
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C         THIS ROUTINE PROVIDES AN APPROXIMATION **RESULT** TO THE
C     INTEGRAL OF AN OVERALL INTEGRAND **HFUN(X)*GPER(X)** OVER A
C     SEMI-INFINITE INTERVAL (AZERO,INFINITY).  THE OVERALL INTEGRAND
C     FUNCTION SHOULD BE ULTIMATELY OSCILLATING IN SIGN AND PERIODIC
C     WITH USER-PROVIDED **PERIOD**.
C
C         THE ROUTINE CONSIDERS A SEQUENCE OF INTERVALS.  ALL THE
C     THE INTERVALS, EXCEPT THE FIRST, ARE OF LENGTH 0.5*PERIOD.
C     IT USES A QUADRATURE RULE TO APPROXIMATE THE INTEGRAL OVER
C     EACH INTERVAL.  THESE APPROXIMATIONS ARE STORED IN
C     QLIST(J)   J = 1,2,3,...   .  THE VALUES OF QLIST ULTIMATELY
C     OSCILLATE IN SIGN.  THEN THE ROUTINE USES A SERIES ACCELERATION
C     TECHNIQUE BASED ON THE EULER TRANSFORMATION TO SUM THIS SEQUENCE.
C
C
C INPUT PARAMETERS:
C
C   AZERO         LOWER LIMIT OF INTEGRATION
C
C   PERIOD        LEAST POSITIVE PERIOD OR ULTIMATE PERIOD OF
C                 INTEGRAND FUNCTION
C
C   RFIRST        RIGHT-HAND ENDPOINT OF FIRST INTERVAL. SUGGESTED VALUE
C                 FOR STRAIGHTFORWARD PROBLEMS IS AZERO.
C                 SEE NOTE 1 BELOW
C
C   EPS           THE REQUESTED ACCURACY
C
C   NQUAD         INTEGER. NQ = ABS(NQUAD) IS THE NUMBER OFABSCISSAS TO
C                 BE USED BY THE QUADRATURE RULE IN EACH INTERVAL.
C                     NQUAD > 1    TRAPEZOIDAL RULE IS USED
C                     NQUAD < 0    RULE SPECIFIED BY GAUSS BELOW IS USED
C
C   NDIM1,NDIM2   DIMENSIONS OF THE OUTPUT ARRAY **WORK**.
C                 NDIM1 = 10 IMPLIES NMAX = 100
C                 NDIM1 .NE. 10 IMPLIES NMAX = MIN(100,NDIM1)
C                 NMAX IS A PHYSICAL LIMIT.  IF CALCULATION IS NOT
C                 COMPLETE AFTER NMAX INTERVALS HAVE BEEN CONSIDERED
C                 IT IS THEN ABANDONED. (SEE NOTE 5 BELOW)
C                 NDIM2  SUGGESTED VALUE IS 15. NORMALLY SHOULD EXCEED 4
C
C   GAUSS         NAME OF A SUBROUTINE WHICH PROVIDES WEIGHTS
C                 AND ABSCISSAS.  WHEN NQUAD > 0 THIS IS NOT CALLED.
C                 SEE NOTE 2 BELOW FOR USE WHEN NQUAD < 0.
C
C   HFUN          NAME OF A FUNCTION SUBROUTINE.  THIS IS ONE FACTOR OF
C                 THE INTEGRAND FUNCTION. THE OTHER FACTOR IS GPER BELOW
C
C   GPER          NAME OF A FUNCTION SUBROUTINE.  THIS IS THE
C                 CO-FACTOR OF HFUN IN THE INTEGRAND FUNCTION.
C                 GPER MUST BE PERIODIC WITH PERIOD COINCIDING WITH THE
C                 INPUT PARAMETER, PERIOD ABOVE. (SEE NOTE 3 BELOW)
C
C
C OUTPUT PARAMETERS:
C
C   WORK(NDIM1,NDIM2)    HOLDS COMPLETED FINITE AVERAGE TABLE
C
C   SAVPER( )         ARRAY TO SAVE FUNCTION VALUES OF GPER.
C                     DIMENSION NOT LESS THAN 2*IABS(NQUAD)
C
C
C   WEIGHT( )         ARRAY TO HOLD WEIGHTS FOR GAUSSIAN QUADRATURE
C                     DIMENSION NOT LESS THAN MAX(1,-NQUAD)
C
C   ABSCIS( )         ARRAY TO HOLD ABSCISSAS FOR GAUSSIAN QUADRATURE
C                     DIMENSION NOT LESS THAN MAX(1,-NQUAD)
C
C   QLIST(100)        ARRAY TO HOLD THE RESULT OF EACH QUADRATURE
C
C   RESULT            OVERALL RESULT OF THE INTEGRATION
C
C   ISTATE(6)         VECTOR OF INTEGERS GIVES STATUS OF RESULT
C
C     ISTATE(1)       AN INDICATOR WARNING ABOUT SUSPICIOUS RESULTS -
C                     ZERO:      THE RUN WAS APPARENTLY SUCCESSFUL.
C                     POSITIVE:  THE RUN WAS APPARENTLY SUCCESSFUL , BUT
C                                THERE ARE UNSATISFACTORY FEATURES OF
C                                POSSIBLE INTEREST TO THE SOPHISTICATED
C                                USER
C                     NEGATIVE:  UNSUCCESSFUL RUN - DISREGARD THE RESULT
C                     SEE NOTES BELOW FOR COMPLETE SPECIFICATION
C     ISTATE(2)       LSIGCH - INDICATES LAST INTERVAL IN WHICH THE SIGN
C                                OF THE INTEGRAL COINCIDED WITH THAT OF
C                                THE INTEGRAL OVER THE PREVIOUS INTERVAL
C                                SEE NOTE 4 BELOW.
C     ISTATE(3)       NOW  - THE NUMBER OF FINITE INTEGRALS, QLIST(Q),
C                                EVALUATED IN THE CALCULATION.
C     ISTATE(4)       NCOL - THE COLUMN OF THE FINITE AVERAGE TABLE,
C                                (WORK), ON WHICH THE  RESULT IS BASED.
C     ISTATE(5)       NROW - THE ROW OF THE FINITE AVERAGE TABLE,
C                                (WORK), ON WHICH THE  RESULT IS BASED.
C     ISTATE(6)       NCOUNT - THE NUMBER OF CALLS TO FUNCTION HFUN.
C
C
C  NOTE 1
C
C  RFIRST
C     THIS ALLOWS THE USER TO LOCATE HIS SUBDIVISION.  FOR CAUTIOUS
C  RUNNING, ARRANGE RFIRST TO COINCIDE WITH AN "ULTIMATE ZERO".  FOR
C  SLIGHTLY MORE ADVENTUROUS BUT LESS RELIABLE RUNNING, ARRANGE
C  RFIRST TO COINCIDE WITH AN "ULTIMATE" PEAK.  OTHERWISE, SET
C  RFIRST < AZERO, IN WHICH CASE, OSCINT USES AZERO INSTEAD OF  RFIRST.
C
C
C  NOTE 2
C
C  GAUSS
C     THIS IS THE NAME OF A USER-PROVIDED SUBROUTINE OF THE FORM
C          GAUSS(ITYPE,A,B,C,D,N,WEIGHT(N),ABS(N)CIS,IFAIL)
C  IT IS CALLED ONLY WHEN NQUAD IS NEGATIVE WITH N = -NQUAD.  IT
C  RETURNS A SET OF WEIGHTS AND ABSCISSAS, SUITABLE FOR INTEGRATION
C  OVER THE INTERVAL (-1,1).
C     THE FIRST FIVE INPUT PARAMETERS OF GAUSS SHOULD BE IGNORED.
C  IFAIL MAY BE USED FOR A WARNING MESSAGE.  IF GAUSS RETURNS
C  IFAIL .NE. 0, OSCINT ABORTS, SETTING ISTATE(1) = -4000.
C     THE USER WHO HAS THE NAG LIBRARY AVAILABLE MAY SET GAUSS = D01BCF.
C
C
C  NOTE 3
C
C  HFUN AND GPER
C     ONE MAY ALWAYS USE HFUN FOR THE INTEGRAND FUNCTION, F(X), AND
C  CODE GPER TO RETURN THE VALUE 1.0D0.  HOWEVER, WHEN F(X) HAS A
C  PERIODIC FACTOR, J(X), REPETITIVE EVALUATION OF J(X) IN EACH
C  INTERVAL MAY BE AVOIDED BY SETTING GPER = J(X) AND HFUN = H(X).
C  THE ROUTINE CHECKS THAT GPER IS INDEED PERIODIC BY MAKING SOME
C  EVALUATIONS IN THE THIRD, FOURTH, AND FIFTH INTERVALS.  IF IT
C  FINDS GPER IS NOT PERIODIC, IT TERMINATES WITH  ISTATE(1) = -5000.
C  IF PERIOD <= 10**-5 , THE ROUTINE TERMINATES WITH ISTATE(1) = -3000.
C
C
C  NOTE 4
C
C  TERMINATION
C     IN THIS NOTE, THE "NORMAL SIGN CHANGE PATTERN" IS ONE IN WHICH
C  SUCCESSIVE VALUES OF QLIST(Q) OSCILLATE IN SIGN.  THE "GRACE
C  PERIOD" IS Q<=10 WHEN THE NORMAL SIGN CHANGE PATTERN IS NOT
C  INSISTED ON.  LSIGCH IS THE HIGHEST VALUE OF Q FOR WHICH
C  QLIST(Q)*QLIST(Q-1) IS POSITIVE.  TERMINATION COMES ABOUT
C  AFTER CALCULATING QLIST(NOW), WHEN EITHER
C
C  (A)  AN APPROXIMATION OF SUFFICIENT ACCURACY IS CURRENTLY AVAIL-
C       ABLE  (THE ROUTINE SETS ISTATE(1) = MAX(0,4-(NOW-LSIGCH))
C  OR
C
C  (B)  THE NORMAL SIGN CHANGE PATTERN IS VIOLATED AFTER THE GRACE
C       PERIOD.   (THE ROUTINE SETS ISTATE(1) = -200)
C   OR
C
C  (C)  THE USER SET LIMIT, NMAX, OF INTERVALS HAVE BEEN CALCULATED
C       I.E.  NOW = NMAX.  (THE ROUTINE SETS ISTATE(1) = -100)
C
C     THE ROUTINE CHECKS (A), THEN (B), AND THEN (C).  TERMINATION
C  UNDER (A) MAY OCCUR WITHOUT ANY NORMAL SIGN CHANGE PATTERN
C  EMERGING.  THIS MAY BE DUE TO MISUSE OF THE ROUTINE, BUT THE
C  PROBLEM IS SMALL ENOUGH TO BE CORRECTLY HANDLED.  IN THIS CASE,
C  ISTATE(1) = MAX(0,4-(NOW-LSIGCH)) MAY BE A SMALL POSITIVE  INTEGER.
C
C
C  NOTE 5
C
C  VARIABLE DIMENSIONS AND STORAGE ECONOMY
C     THE PROGRAM WHICH CALLS OSCINT HAS TO PROVIDE NUMERICAL VALUES
C  OF THE INPUT PARAMETERS NQUAD, NDIM1 AND NDIM2. IT MUST ALSO
C  INCLUDE A DIMENSION STATEMENT IN WHICH THE FIRST FIVE AND THE
C  LAST OUTPUT PARAMETERS OF OSCINT ARE DIMENSIONED.  NOTE THAT THE
C  DIMENSIONS OF WEIGHT AND ABSCIS ARE BOTH AT LEAST 1 WHEN NQUAD
C  IS POSITIVE AND AT LEAST -NQUAD OTHERWISE.  HOWEVER, THE
C  DIMENSION OF SAVPER IS AT LEAST 2*ABS(NQUAD) REGARDLESS OF THE
C  SIGN OF NQUAD.
C     GENERALLY, THESE VARIABLE DIMENSION STATEMENTS ALLOW ECONOMIC
C  USE OF STORAGE.  THERE IS ONE FURTHER STORAGE SAVING FEATURE.
C  WHEN NDIM1 = 10, THE ROUTINE CARRIES OUT THE SAME CALCULATION
C  AS IT WOULD IF NDIM1 = 100.  HOWEVER, INSTEAD OF USING A WORK
C  ARRAY OF DIMENSION (100,NDIM2), IT USES A WORK ARRAY OF DIMENSION
C  (10,NDIM2) AND OVERWRITES IT, AS AND WHEN NECESSARY, AS THE
C  CALCULATION PROCEEDS. IN THIS CASE, OSCINT BEHAVES AS IF NDIM2 WERE
C  REPLACED BY  MIN(20,NDIM2).
C     SOME OBVIOUS ERRORS IN DIMENSIONING, SUCH AS NEGATIVE
C  DIMENSIONS, CAUSE THE ROUTINE TO TERMINATE WITH
C         ISTATE(1) = -6000
C  HOWEVER, INADEQUATE DIMENSIONING IN THE CALLING PROGRAM MAY GO
C  UNDETECTED AND MAY LEAD TO RANDOM OR CHAOTIC OUTPUT.
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C      IMPLICIT NONE
C     *** T. WIEDER, 15.04.1998 ***
C      DOUBLE PRECISION SAVPER(*),WEIGHT(*),ABSCIS(*)
C
C     .. Parameters ..
      INTEGER NDIMMA
      PARAMETER (NDIMMA=1000)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION AZERO,EPS,PERIOD,RESULT,RFIRST
      INTEGER NDIM1,NDIM2,NQUAD
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION ABSCIS(36),QLIST(0:NDIMMA),SAVPER(0:256),
     +                 WEIGHT(36),WORK(0:NDIM1,0:NDIM2)
      INTEGER ISTATE(6)
C     ..
C     .. Function Arguments ..
      DOUBLE PRECISION GPER,HFUN
      EXTERNAL GPER,HFUN
C     ..
C     .. Subroutine Arguments ..
      EXTERNAL GAUSS
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION CURR,HASPER,PREV,S,WMIN
      INTEGER I,II,ISUB,J,JP,K,MJ,MJ1,MJ2,MJK,MJM2,NMAX,NML,NOW,NOWJP,
     +        NROUND,NROW
C     ..
C     .. External Functions ..
      DOUBLE PRECISION QRULE
      EXTERNAL QRULE
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,DABS,MAX,MIN,MOD
C     ..
      NROUND = 0
      HASPER = 0.5D0*PERIOD
      WMIN = 1.0D0
      DO 10 I = 1,6
          ISTATE(I) = 0
   10 CONTINUE
      IF (NQUAD.EQ.0 .OR. NQUAD.EQ.1 .OR. NDIM1.LT.1 .OR.
     +    NDIM2.LT.1) ISTATE(1) = -6000
      IF (PERIOD.LT.10.0D-5) ISTATE(1) = -5000
      IF (ISTATE(1).LT.0) GO TO 120
C
C     *** T. WIEDER, 18.02.1998 ***
      IF (NDIM1.EQ.10) THEN
          NMAX = 100

      ELSE
C          NMAX = MIN(100,NDIM1)
          NMAX = MIN(NDIMMA,NDIM1)
      END IF
C
      JP = 0
      S = 0.0D0
C
C      LOOP TO CONSTRUCT TABLE
C
   20 IF (JP.EQ.NMAX) ISTATE(1) = -100
      IF (ISTATE(1).NE.0) GO TO 40
C
      K = 0
      JP = JP + 1
      J = MIN(JP,NDIM2)
      IF (NDIM1.EQ.10) J = MIN(J,20)
   30 IF (K.LT.J) THEN
          K = K + 1
          IF (NDIM1.EQ.10) THEN
              MJ = MOD(JP,10)
              MJM2 = MOD(JP-2,10)
              MJ1 = MOD(JP-K+1,10)
              MJ2 = MOD(JP-K+2,10)
              MJK = MOD(JP-K,10)
              IF (MJ.EQ.0) MJ = 10
              IF (MJM2.EQ.0) MJM2 = 10
              IF (MJ1.EQ.0) MJ1 = 10
              IF (MJ2.EQ.0) MJ2 = 10
              IF (MJK.EQ.0) MJK = 10
C
          ELSE
              MJ = JP
              MJM2 = JP - 2
              MJ1 = JP - K + 1
              MJ2 = JP - K + 2
              MJK = JP - K
          END IF
C
          IF (K.EQ.1) THEN
C
C      *** AMENDMENT SUGGESTED BY DR. T.R. HOPKINS, 15.04.1998 ***
C      WORK(MJ,K)=QRULE(JP-1,AZERO,HASPER,RFIRST,NQUAD,NDIM1,NDIM2,GAUSS
C     1HFUN,GPER,SAVPER,WEIGHT,ABSCIS,QLIST,ISTATE)
              WORK(MJ,K) = QRULE(JP-1,AZERO,HASPER,RFIRST,NQUAD,GAUSS,
     +                     HFUN,GPER,SAVPER,WEIGHT,ABSCIS,QLIST,ISTATE)
C
          ELSE
              WORK(MJ1,K) = (WORK(MJ1,K-1)+WORK(MJ2,K-1))/2.0D0
              IF (DABS(WORK(MJ1,K)).LT.WMIN) THEN
                  WMIN = DABS(WORK(MJ1,K))
                  NOW = K
                  NROW = MJ1
                  NOWJP = JP
              END IF
C
              CURR = ABS(WORK(MJ1,K))
              PREV = ABS(WORK(MJK,K))
              IF (JP.NE.K) THEN
                  IF ((CURR.LT.EPS) .AND. (PREV.LT.EPS)) THEN
                      NOW = K
                      NROW = MJ1
                      NOWJP = JP
                      GO TO 40
C
                  END IF
C
              END IF
C
          END IF
C
          GO TO 30
C
      END IF
C
      GO TO 20
C
   40 DO 50 I = 1,NOWJP - NOW
          S = S + QLIST(I)
   50 CONTINUE
      IF (NDIM1.NE.10) GO TO 100
      NROUND = NOW - 10
      IF (NROUND.LE.0) NROUND = 0
      IF (NROUND.LE.0) GO TO 100
C
      DO 60 I = NOWJP - NOW + 1,NOWJP - 10
          S = S + QLIST(I)
   60 CONTINUE
      ISUB = NOWJP - 9
   70 IF (ISUB.GT.10) ISUB = ISUB - 10
      IF (ISUB.GT.10) GO TO 70
      DO 80 J = 1,NROUND
          S = S + WORK(ISUB,J)/2.0D0
   80 CONTINUE
      DO 90 I = NROW,NROW + NROUND - 1
          II = I
          IF (I.GT.10) II = I - 10
          S = S - WORK(II,NROUND+1)
   90 CONTINUE
C
  100 DO 110 J = NROUND + 1,NOW - 1
          S = S + WORK(NROW,J)/2.0D0
  110 CONTINUE
      S = S + WORK(NROW,NOW)
      RESULT = S
      ISTATE(3) = JP
      ISTATE(4) = NOW
      ISTATE(5) = NROW
      NML = ISTATE(3) - ISTATE(2)
      IF (NML.LT.4 .AND. ISTATE(1).EQ.0) ISTATE(1) = MAX(0,4-NML)
  120 RETURN

      END
C     END OF SUBROUTINE *** OSCINT ***
C
C     BEGIN OF FUNCTION *** QRULE ***
C
C      *** AMENDMENT SUGGESTED BY DR. T.R. HOPKINS, 15.04.1998 ***
C      DOUBLE PRECISION FUNCTION QRULE (J,AZERO,HASPER,RFIRST,NQUAD,
C     1NDIM1,NDIM2,GAUSS,HFUN,GPER,SAVPER,WEIGHT,ABSCIS,QLIST,ISTATE)
      DOUBLE PRECISION FUNCTION QRULE(J,AZERO,HASPER,RFIRST,NQUAD,GAUSS,
     +                                HFUN,GPER,SAVPER,WEIGHT,ABSCIS,
     +                                QLIST,ISTATE)
C
C     ##################################################################
C
C     *** T. WIEDER, 05.03.1998 ***
C
C     SOURCE FOR *** FUNCTION QRULE ***:
C
C     SUBPROGRAM LIBRARY TOMS
C
C     AVAILABILITIY: PUBLIC DOMAIN
C
C     DEVELOPER: JAMES LYNESS AND GWENDOLEN HINES.
C
C     DISTRIBUTOR: TRANSACTIONS ON MATHEMATICAL SOFTWARE,
C     VOL. 12 (1986), PP. 24 - 25.
C
C     ##################################################################
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C          THIS ROUTINE EVALUATES THE INTEGRAL OF THE FUNCTION,
C     HFUN(X)*GPER(X) OVER THE INTERVAL (A,B) WHERE:
C
C          GENERALLY (J>0)     A = PREVIOUS VALUE OF B
C                              B = A + HASPER
C                    (J=0)     A = AZERO
C                              B = RFIRST WHEN RFIRST > AZERO
C                              B = AZERO + HASPER WHEN RFIRST <= AZERO
C
C
C      INPUT PARAMETERS:
C
C               J                  DEFINES WHICH TERM, QLIST(J), OF
C                                  THE SERIES IS BEING EVALUATED
C
C               HASPER             HALF THE PERIOD
C
C
C      OTHER INPUT AND OUTPUT PARAMETERS ARE IDENTICAL TO THOSE IN
C      OSCINT (DESCRIBED ABOVE).
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     *** T. WIEDER, 18.02.1998 ***
C      IMPLICIT NONE
C     ORIGINAL VALUE FOR *** IGRACE *** = 9
C     ORIGINAL VALUE FOR *** NDIMMA *** = 100
C     *** T. WIEDER, 15.04.1998 ***
C      DOUBLE PRECISION SAVPER(*),WEIGHT(*),ABSCIS(*)
C
C     .. Parameters ..
      INTEGER IGRACE,NDIMMA
      PARAMETER (IGRACE=99,NDIMMA=1000)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION AZERO,HASPER,RFIRST
      INTEGER J,NQUAD
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION ABSCIS(36),QLIST(0:NDIMMA),SAVPER(0:256),
     +                 WEIGHT(36)
      INTEGER ISTATE(6)
C     ..
C     .. Function Arguments ..
      DOUBLE PRECISION GPER,HFUN
      EXTERNAL GPER,HFUN
C     ..
C     .. Subroutine Arguments ..
      EXTERNAL GAUSS
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION A,AA,B,BB,CC,DD,DIFF,FUN,GMAX,TENDPT,WSUM,WT,XI,Y
      INTEGER I,IELM,IFAIL,INDEX,ITYPE,NPTS
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,DABS,MOD
C     ..
C     .. Save statement ..
      SAVE B,TENDPT
C     ..
      NPTS = ABS(NQUAD)
C
      IF (J.EQ.0) THEN
          INDEX = 0
          DO 10 I = 1,NDIMMA
              QLIST(I) = 0.0D0
   10     CONTINUE
C        CALL GAUSS IF NQUAD < 0
          IF (NQUAD.LT.0) THEN
              ITYPE = 0
              AA = -1.0D0
              BB = 1.0D0
              CC = 0.0D0
              DD = 0.0D0
              IFAIL = 0
              CALL GAUSS(ITYPE,AA,BB,CC,DD,NPTS,WEIGHT,ABSCIS,IFAIL)
              IF (IFAIL.NE.0) THEN
                  ISTATE(1) = -4000
                  GO TO 40
C
              END IF
C
          END IF
C
      END IF
C
      QRULE = 0.0D0
C
C     SET INTERVAL ENDPOINTS, A AND B
      IF (J.NE.0) THEN
          A = B
C
      ELSE
          A = AZERO
      END IF
C
      IF (RFIRST.LE.AZERO .OR. J.NE.0) THEN
          B = A + HASPER
C
      ELSE
          B = RFIRST
      END IF
C
C     LOOP ON ABSCISSAS FOR INTERVAL #J - STARTS TO CALCULATE QRULE
      DO 30 I = 1,NPTS
          XI = I
C
C        CALCULATE ABSCISSA Y
          IF (NQUAD.LT.0) THEN
              Y = (B-A)*ABSCIS(I)/2.0D0 + (B+A)/2.0D0
              WT = WEIGHT(I)
C
          ELSE
              Y = ((XI-1)*B+A* (NPTS-I))/ (NPTS-1)
              WT = 1.0D0
          END IF
C
C        IELM IS LOCATION IN SAVPER FOR GPER FUNCTION VALUES
C        J=0: IGNORES. J=1,2: WHERE TO PUT VALUE. J>2: WHERE TO GET
C        VALUE FROM
          IELM = NPTS*MOD(J-1,2) + I
C
          IF (J.GT.0 .AND. I.EQ.1 .AND. NQUAD.GT.0) GO TO 20
          IF (J.EQ.0) THEN
              FUN = HFUN(Y)*GPER(Y)
C
          ELSE IF (J.LT.3) THEN
              SAVPER(IELM) = GPER(Y)
C           CHECK FOR CONSTANT FUNCTION
              IF (SAVPER(IELM).EQ.SAVPER(1)) INDEX = INDEX + 1
C     CODE MODIFIED BY DR. T.R. HOPKINS, 28.04.1998:
              IF (IELM.EQ.1) THEN
                  GMAX = DABS(SAVPER(1))

              ELSE
                  IF (DABS(SAVPER(IELM)).GT.DABS(SAVPER(IELM-1))) THEN
                      GMAX = DABS(SAVPER(IELM))
                  END IF

              END IF
C
              FUN = HFUN(Y)*SAVPER(IELM)
C
          ELSE IF (J.LT.5) THEN
C           CHECK THAT GPER IS PERIODIC WITH PERIOD=PERIOD
              DIFF = DABS(SAVPER(IELM)-GPER(Y))
C
C     *** T. WIEDER, 19.02.1998 ***
C      IF (DIFF.LT.GMAX*1.0D-5) THEN
              IF (DIFF.LT.GMAX*1.0D1) THEN
                  FUN = HFUN(Y)*SAVPER(IELM)
C
              ELSE
                  ISTATE(1) = -3000
              END IF
C
          ELSE
              IF (INDEX.EQ.2*NPTS) THEN
                  FUN = HFUN(Y)
C
              ELSE
                  FUN = HFUN(Y)*SAVPER(IELM)
              END IF
C
          END IF
C
          ISTATE(6) = ISTATE(6) + 1
C
   20     IF (NQUAD.GT.0) THEN
              IF (I.EQ.1) THEN
                  IF (J.EQ.0) FUN = FUN/2.0D0
                  IF (J.GT.0) FUN = TENDPT
              END IF
C
              IF (I.EQ.NPTS) THEN
                  FUN = FUN/2.0D0
                  TENDPT = FUN
              END IF
C
              QRULE = QRULE + FUN
C
          ELSE
              QRULE = QRULE + WT*FUN
          END IF
C
   30 CONTINUE
C     LOOP ON ABSCISSA FOR INTERVAL #J ENDS
      IF (INDEX.EQ.2*NPTS .AND. SAVPER(1).NE.1) THEN
          QRULE = QRULE*SAVPER(1)
      END IF
C
      IF (NQUAD.GT.0) WSUM = NPTS - 1
      IF (NQUAD.LT.0) WSUM = 2.0D0
      QRULE = QRULE* (B-A)/WSUM
      QLIST(J+1) = QRULE
C
      IF (QRULE*QLIST(J).GT.0) ISTATE(2) = J
C
C     *** T. WIEDER, 18.02.1998 ***
C      IF (ISTATE(2).GT.9) ISTATE(1) = -200
      IF (ISTATE(2).GT.IGRACE) ISTATE(1) = -200
C
   40 RETURN

      END
C     END OF SUBROUTINE *** QRULE ***
C
C     BEGIN OF SUBROUTINE *** G5AND9 ***
      SUBROUTINE G5AND9(ITYPE,A,B,C,D,NQUAD,WEIGHT,ABSCIS,IERR)
C
C     ##################################################################
C
C     *** T. WIEDER, 05.03.1998 ***
C
C     SOURCE FOR SUBROUTINE *** G5AND9 ***:
C
C     SUBPROGRAM LIBRARY TOMS
C
C     AVAILABILITIY: PUBLIC DOMAIN
C
C     DEVELOPER: JAMES LYNESS AND GWENDOLEN HINES.
C
C     DISTRIBUTOR: TRANSACTIONS ON MATHEMATICAL SOFTWARE,
C     VOL. 12 (1986), PP. 24 - 25.
C
C     ##################################################################
C
C      THIS IS AN EXTRACT FROM A QUADRATURE ROUTINE CONSTRUCTED ONLY FOR
C      USE IN A DRIVER WHICH ILLUSTRATES OSCINT. A NAG LIBRARY
C      SUBSCRIBER MAY REPLACE THIS BY D01BCF FOR GENERAL USE.
C
C      IMPLICIT NONE
C
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION A,B,C,D
      INTEGER IERR,ITYPE,NQUAD
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION ABSCIS(NQUAD),WEIGHT(NQUAD)
C     ..
C     .. Local Scalars ..
      INTEGER I
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION ABSC5(5),ABSC9(9),WT5(5),WT9(9)
C     ..
C     .. Data statements ..
C      THE FOLLOWING DATA ARE WEIGHTS AND ABSCISSAS OF THE FIVE POINT
C      AND THE NINE POINT GAUSS-LEGENDRE QUADRATURE RULES RESPECTIVELY.
C      NORMALISED TO THE INTERVAL (-1,1).
      DATA WT5/.2369268850561891D0,.4786286704993665D0,
     +     .5688888888888889D0,.4786286704993665D0,.2369268850561891D0/
      DATA ABSC5/-.9061798459386640D0,-.5384693101056831D0,0.0D0,
     +     .5384693101056830D0,.9061798459386640D0/
      DATA WT9/.0812743883615745D0,.1806481606948574D0,
     +     .2606106964029355D0,.3123470770400029D0,.3302393550012598D0,
     +     .3123470770400029D0,.2606106964029356D0,.1806481606948575D0,
     +     .8127438836157467D-01/
      DATA ABSC9/-.9681602395076261D0,-.8360311073266358D0,
     +     -.6133714327005904D0,-.3242534234038090D0,0.0D0,
     +     .3242534234038087D0,.6133714327005902D0,.8360311073266357D0,
     +     .9681602395076260D0/
C     ..
C
C     *** T. WIEDER, 16.04.1998 ***
C     THIS IS A DIRTY TRICK TO PREVENT COMPILER FROM WARNING MESSAGE:
      A = A
      B = B
      C = C
      D = D
      ITYPE = ITYPE
C
      IERR = 59
C
      IF (NQUAD.EQ.5) THEN
          IERR = 0
          DO 10 I = 1,5
              ABSCIS(I) = ABSC5(I)
              WEIGHT(I) = WT5(I)
   10     CONTINUE
      END IF
C
      IF (NQUAD.EQ.9) THEN
          IERR = 0
          DO 20 I = 1,9
              ABSCIS(I) = ABSC9(I)
              WEIGHT(I) = WT9(I)
   20     CONTINUE
      END IF
C
      RETURN

      END
C     END OF SUBROUTINE *** G5AND9 ***
C     BEGIN OF SUBROUTINE *** TBESSE ***
      SUBROUTINE TBESSE
C
C     #################################################################
C
C     FIRST VERSION:    06.02.1998
C     PREVIOUS VERSION: 06.02.1998
C     LATEST VERSION:   05.03.1998
C
C     PROGRAMMED BY:
C     THOMAS WIEDER
C     E-MAIL: WIEDER@HRZ.UNI-KASSEL.DE
C
C     PURPOSE:
C
C     TEST SUBROUTINE *** RJBESL ***.
C
C     #################################################################
C
C      IMPLICIT NONE
C
C     .. Local Scalars ..
      DOUBLE PRECISION ALPHA,NU,X
      INTEGER N,NB,NCALC
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION RESULT(101)
C     ..
C     .. External Subroutines ..
      EXTERNAL RJBESL
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC INT
      integer i1mach, nout, nin
      nout = i1mach(2)
      nin = i1mach(1)
C     ..
      OPEN (90,FILE='TBESSE.out',STATUS='UNKNOWN',ACCESS='SEQUENTIAL',
     +     ERR=30)
      WRITE(NOUT,FMT=9000)
      READ(NIN,FMT=*,ERR=20) X,NU,ALPHA
C
      NB = INT(NU) + 1
      WRITE (90,FMT=9010)
      WRITE(NOUT,FMT=9010)
C
      CALL RJBESL(X,ALPHA,NB,RESULT,NCALC)
C
      WRITE(NOUT,FMT=9020) NCALC
      WRITE(NOUT,FMT=9030) NB
      IF (NCALC.NE.NB) THEN
          WRITE(NOUT,FMT=9040)
      END IF

      DO 10 N = 1,NB
          WRITE (90,FMT=9050) N - 1 + ALPHA,X,RESULT(N)
          WRITE(NOUT,FMT=9050) N - 1 + ALPHA,X,RESULT(N)
   10 CONTINUE
C
      CLOSE (90)
      GO TO 40

   20 WRITE(NOUT,FMT=9060)
      GO TO 40

   30 WRITE(NOUT,FMT=9070)
   40 STOP

 9000 FORMAT (' MESSAGE TBESSE:',/,' GIVE X, NU, ALPHA',:,/)
 9010 FORMAT ('# MESSAGE TBESSE: NU+ALPHA, X, BESSEL(NU+ALPHA,X)')
 9020 FORMAT ('# MESSAGE TBESSE: NCALC =',I6)
 9030 FORMAT ('# MESSAGE TBESSE: NB =',I6)
 9040 FORMAT ('# MESSAGE TBESSE: NCALC NOT EQUAL NB!')
 9050 FORMAT (D17.9,2X,D17.9,2X,D17.9)
 9060 FORMAT (/,' MESSAGE TBESSE: ERROR DURING INPUT!')
 9070 FORMAT (/,'MESSAGE TBESSE: ERROR ON OPENING FILE!')
      END
C     END OF SUBROUTINE *** TBESSE ***
C
C     BEGIN OF SUBROUTINE *** TINTHP ***
      SUBROUTINE TINTHP
C
C     #################################################################
C
C     FIRST VERSION: 05.02.1998
C     PREVIOUS VERSION: 05.02.1998
C     LATEST VERSION: 06.02.1998
C
C     PROGRAMMED BY:
C     THOMAS WIEDER
C     E-MAIL: WIEDER@HRZ.UNI-KASSEL.DE
C
C     PURPOSE:
C     TEST SUBROUTINE *** INTHP ***.
C
C     #################################################################
C
C      IMPLICIT NONE
C
C     .. Scalars in Common ..
      DOUBLE PRECISION EXACT,XSPLIT
      INTEGER IINT
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION A,B,DCLASS,DX,EPS,PCLASS,RESULT,X
      INTEGER I,INF,MCALL
C     ..
C     .. External Functions ..
      DOUBLE PRECISION TFUNC
      EXTERNAL TFUNC
C     ..
C     .. External Subroutines ..
      EXTERNAL INTHP
C     ..
C     .. Common blocks ..
      COMMON /LOSUNG/XSPLIT,EXACT,IINT
      integer i1mach, nout, nin
      nout = i1mach(2)
      nin = i1mach(1)
C     ..
      OPEN (81,FILE='test.plo',STATUS='unknown',ACCESS='sequential')
C
   10 WRITE(NOUT,FMT=9000)
      READ(NIN,FMT=*,ERR=10) IINT
      IF (IINT.EQ.0) GO TO 30
C
      X = 0.0D0
      DX = 1.0D0
      DO 20 I = 1,1000
          WRITE (81,FMT=*) X,TFUNC(X)
          X = X + DX
   20 CONTINUE
C
C     DO AN INDEFINITE INTEGRAL FROM *** BOUND *** TO INFINITY:
      A = 0.0D0
      B = 1.0D44
      MCALL = 20000
C
C     *** 05.02.1998 ***
C     GOOD VALUE FOR *** EPS *** = 1.0D-3
      EPS = 1.0D-03
C     *** 05.02.1998 ***
C
C     *** 05.02.1998 ***
C     GOOD VALUE FOR *** DCLASS *** = 10.0D0
C     ALL OTHER VALUES RESULT IN FAILORE!
      DCLASS = 10.0D0
C     *** 05.02.1998 ***
C
C     *** 05.02.1998 ***
C     GOOD VALUE FOR *** PCLASS *** = 0.0D0
C     ALL OTHER VALUES RESULT IN FAILORE!
      PCLASS = 0.0D0
C     *** 05.02.1998 ***
C
C     *** 05.02.1998 ***
C     TAKE *** INF = 2 *** FOR NON-OSCILLATING INTEGRALS!
C     TAKE *** INF = 3 *** FOR OSCILLATING INTEGRALS!
      INF = 3
C     *** 05.02.1998 ***
C
      CALL INTHP(A,B,DCLASS,TFUNC,MCALL,PCLASS,EPS,INF,RESULT)
C
      WRITE(NOUT,FMT=*) 'RESULT =',RESULT
      WRITE(NOUT,FMT=*) 'NUMBER OF FUNCTION CALLS MACLL =',MCALL
      WRITE(NOUT,FMT=*) 'INFORMATION FROM INTHP =',INF
      WRITE(NOUT,FMT=*) 'EXACT SOLUTION:',EXACT
C
      GO TO 10
C
   30 RETURN

 9000 FORMAT (1X,' GIVE IDENTIFICATION NUMBER OF INTEGRAND!',/,
     +       ' (1,2,3 4,5,6,7,8,9, 0 = STOP):',/)
      END
