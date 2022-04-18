c
c
      PROGRAM DTEST6
c
c This program is a double-precision version of  test6
c
c      parameter(ak=1.5d0)
c      parameter(ak=2.5d0)
C     .. Parameters ..
      INTEGER NMAX,NCAPM
      PARAMETER (NMAX=10,NCAPM=200)
      DOUBLE PRECISION THETA
      PARAMETER (THETA=1.D-4)
      DOUBLE PRECISION ETA
      PARAMETER (ETA=-1.D0)
      DOUBLE PRECISION AK
      PARAMETER (AK=.5D0)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION CONST,DEN,EPS,EXACT,PI,SUM
      REAL ERR
      INTEGER IERAB,IERGQ,IERR,IROUT,K,KOUNT,M,MU,N,NCAP
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION A(NCAPM),ALPHA(NMAX+1),B(NCAPM),BE(NMAX+1),
     +                 BETA(NMAX+1),E(NCAPM),P0(NCAPM),P1(NCAPM),
     +                 P2(NCAPM),W(NCAPM),WG(NMAX),X(NCAPM),XII(2*NMAX),
     +                 XIR(2*NMAX),ZG(NMAX)
      INTEGER IS(2*NMAX)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION D1MACH,F
      EXTERNAL D1MACH,F
C     ..
C     .. External Subroutines ..
      EXTERNAL DABMOD,DGQRAT,DRECUR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,ATAN,DBLE,REAL
      integer i1mach, nout
      nout = i1mach(2)
C     ..
      IROUT = 1
c      irout=0
      EPS = D1MACH(3)*1.D2
      CALL DRECUR(NCAPM,7,AK-1.D0,0.D0,A,B,IERR)
      PI = 4.D0*ATAN(1.D0)
      WRITE(NOUT,FMT=9000)
      DO 30 N = 1,NMAX
          M = 2*N - 1
c        m=2*((n+1)/2)-1
c        m=1
c        m=0
          IF (M.GE.3) THEN
              DO 10 MU = 1,M - 2,2
                  DEN = ETA**2 + DBLE((MU+1)**2)* (PI**2)
                  XIR(MU) = -ETA/DEN
                  XIR(MU+1) = XIR(MU)
                  XII(MU) = DBLE(MU+1)*PI/DEN
                  XII(MU+1) = -XII(MU)
                  IS(MU) = 1
                  IS(MU+1) = 1
   10         CONTINUE
          END IF

          IF (M.GT.0) THEN
              XIR(M) = -1.D0/ETA
              XII(M) = 0.D0
              IS(M) = 1
          END IF

          CALL DABMOD(N+1,NCAPM,M,EPS,IROUT,A,B,XIR,XII,IS,ALPHA,BETA,
     +                NCAP,KOUNT,IERAB,BE,X,W,E,P0,P1,P2)
          IF (IERAB.NE.0) WRITE(NOUT,FMT=9010) IERAB
          CALL DGQRAT(N,M,ALPHA,BETA,XIR,XII,IS,ZG,WG,CONST,IERGQ,E)
          IF (IERGQ.NE.0) WRITE(NOUT,FMT=9020) IERGQ
          SUM = 0.D0
          DO 20 K = 1,N
              SUM = SUM + WG(K)*F(ZG(K),ETA,THETA)
   20     CONTINUE
          EXACT = .379708865998074D0
c        exact=.52608888707965d0
c        exact=1.266569126543118d0
          ERR = ABS(REAL((SUM-EXACT)/EXACT))
          WRITE(NOUT,FMT=9030) N,M,SUM,ERR,CONST
   30 CONTINUE
      STOP

 9000 FORMAT (5X,'n',4X,'m',8X,'integral',10X,'rel error',4X,
     +       'err const',/)
 9010 FORMAT (1X,'ierab=',I3)
 9020 FORMAT (1X,'iergq=',I3)
 9030 FORMAT (1X,2I5,D23.14,E13.4,D13.4)
      END

      DOUBLE PRECISION FUNCTION F(T,ETA,THETA)
C     .. Scalar Arguments ..
      DOUBLE PRECISION ETA,T,THETA
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION S,S1,TERM,X
      INTEGER L
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,DBLE,EXP,SQRT
      integer i1mach, nout
C     ..
      X = T - ETA
      IF (ABS(X).LE.1.D0) THEN
          L = 0
          TERM = X
          S1 = TERM
   10     S = S1
          L = L + 1
          IF (L.GT.200) GO TO 20
          TERM = X*TERM/DBLE(L+1)
          S1 = S + TERM
          IF (S1.NE.S) GO TO 10
          F = T*SQRT(1.D0+.5D0*THETA*T)/ (EXP(-T)*S)
          RETURN

      ELSE
          F = T*SQRT(1.D0+.5D0*THETA*T)/ (EXP(-ETA)-EXP(-T))
          RETURN

      END IF

   20 nout = i1mach(2)
      WRITE(NOUT,FMT=9000)
      RETURN

 9000 FORMAT (1X,'exp series does not converge')
      END
