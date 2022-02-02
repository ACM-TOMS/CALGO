c
c
      PROGRAM TEST7
c
c This program does the same as  test6  but treats the real pole	
c close to the origin separately.
c
c
c Be sure that ncapm is greater or equal to kapmax
c
c      parameter(ak=1.5)
c      parameter(ak=2.5)
C     .. Parameters ..
      INTEGER NMAX,NCAPM,KAPMAX
      REAL THETA
      PARAMETER (NMAX=10,NCAPM=50,KAPMAX=20,THETA=1.e-4)
      REAL AK
      PARAMETER (AK=.5)
C     ..
C     .. Local Scalars ..
      COMPLEX Z
      REAL C1,C2,CONST,EPS,ERR,ETA,EXACT,FMU,FMU2,FNU,FNU2,HN,HP,P,PI,
     +     SQPI,SUM,SUM0,T,X1,XA
      INTEGER IERAB,IERGC,IERGQ,IERK,IERR,IROUT,J,K,KOUNT,M,MU,N,NCAP,
     +        NU,NU0,NU1
C     ..
C     .. Local Arrays ..
      COMPLEX RHO(1),ROLD(1)
      REAL A(NCAPM),A1(NMAX+1),ALPHA(NMAX+1),B(NCAPM),B1(NMAX+1),
     +     BE(NMAX+1),BETA(NMAX+1),E(NCAPM),P0(NCAPM),P1(NCAPM),
     +     P2(NCAPM),W(NCAPM),WG(NMAX),X(NCAPM),XII(2*NMAX),XIR(2*NMAX),
     +     ZG(NMAX)
      INTEGER IS(2*NMAX)
C     ..
C     .. External Functions ..
      REAL F,R1MACH
      EXTERNAL F,R1MACH
C     ..
C     .. External Subroutines ..
      EXTERNAL ABMOD,GCHRIS,GQRAT,KNUM,RECUR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,ATAN,CMPLX,EXP,REAL,SQRT
      integer i1mach, nout
      nout = i1mach(2)
C     ..
      PI = 4.*ATAN(1.)
      SQPI = SQRT(PI)
      C1 = PI
c      c1=-pi
c      c1=pi
      C2 = SQPI
c      c2=.5*sqpi
c      c2=.75*sqpi
      IROUT = 1
c      irout=0
      EPS = R1MACH(3)*1.e2
      CALL RECUR(NCAPM,7,AK-1.,0.,A,B,IERR)
      WRITE(NOUT,FMT=9000)
      DO 60 N = 1,NMAX
          M = 2*N - 1
c        m=2*((n+1)/2)-1
c        m=1
          IF (M.GE.3) THEN
              DO 10 MU = 1,M - 2,2
                  XIR(MU) = 0.
                  XIR(MU+1) = 0.
                  XII(MU) = 1./ (REAL(MU+1)*PI)
                  XII(MU+1) = -XII(MU)
                  IS(MU) = 1
                  IS(MU+1) = 1
   10         CONTINUE
          END IF

          CALL ABMOD(N+1,NCAPM,M-1,EPS,IROUT,A,B,XIR,XII,IS,A1,B1,NCAP,
     +               KOUNT,IERAB,BE,X,W,E,P0,P1,P2)
          IF (IERAB.NE.0) WRITE(NOUT,FMT=9010) IERAB
          ETA = -.001
          XIR(M) = -1./ETA
          XII(M) = 0.
          IS(M) = 1
          X1 = ETA
          XA = ABS(ETA)
          T = 1.
          SUM0 = 0.
          SUM = -1./ (AK-1.)
          J = 0
   20     J = J + 1
          IF (J.GT.100) THEN
              WRITE(NOUT,FMT=9020)
              GO TO 60

          END IF

          T = -XA*T/REAL(J)
          SUM = SUM + T/ (REAL(J+1)-AK)
          IF (SUM.NE.SUM0) THEN
              SUM0 = SUM
              GO TO 20

          END IF

          HP = -EXP(XA)* (C1* (XA** (AK-1.))-C2*SUM)
          HN = 0.
          IF (M.GE.3) THEN
              SUM = 0.
              DO 40 NU = 1, (M-1)/2
                  FNU = REAL(NU)
                  FNU2 = FNU**2
                  P = 1.
                  DO 30 MU = 1, (M-1)/2
                      IF (MU.EQ.NU) GO TO 30
                      FMU = REAL(MU)
                      FMU2 = FMU**2
                      P = (FMU2/ (FMU2-FNU2))*P
   30             CONTINUE
                  NU0 = 10
                  Z = CMPLX(0.,2.*FNU*PI)
                  CALL KNUM(0,NU0,KAPMAX,Z,EPS,A,B,RHO,NU1,IERK,ROLD)
                  IF (IERK.NE.0) WRITE(NOUT,FMT=9030) IERK
                  SUM = SUM + 2.*FNU*PI*P* (2.*FNU*PI*HP-
     +                  REAL(CMPLX(2.*FNU*PI,-X1)*RHO(1)))/
     +                  (X1**2+4.*FNU2* (PI**2))
   40         CONTINUE
              HP = SUM
          END IF

          CALL GCHRIS(N+1,1,A1,B1,X1,HP,HN,ALPHA,BETA,IERGC)
          IF (IERGC.NE.0) WRITE(NOUT,FMT=9040) IERGC
          BETA(1) = -ETA*BETA(1)
          CALL GQRAT(N,M,ALPHA,BETA,XIR,XII,IS,ZG,WG,CONST,IERGQ,E)
          IF (IERGQ.NE.0) WRITE(NOUT,FMT=9050) IERGQ
          SUM = 0.
          DO 50 K = 1,N
              SUM = SUM + WG(K)*F(ZG(K),ETA,THETA)
   50     CONTINUE
          EXACT = 2.2171501
c        exact=1.7800123
c        exact=3.7403845
          ERR = ABS((SUM-EXACT)/EXACT)
          WRITE(NOUT,FMT=9060) N,M,SUM,ERR,CONST
   60 CONTINUE
      STOP

 9000 FORMAT (5X,'n',4X,'m',6X,'integral',7X,'rel error',5X,'err const',
     +       /)
 9010 FORMAT (1X,'ierab=',I3)
 9020 FORMAT (1X,'power series for gamma does not converge')
 9030 FORMAT (1X,'ierk for complex z =',I3)
 9040 FORMAT (1X,'iergc=',I3)
 9050 FORMAT (1X,'iergq=',I3)
 9060 FORMAT (1X,2I5,E17.7,2E14.4)
      END

      REAL FUNCTION F(T,ETA,THETA)
C     .. Scalar Arguments ..
      REAL ETA,T,THETA
C     ..
C     .. Local Scalars ..
      REAL S,S1,TERM,X
      INTEGER L
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,EXP,REAL,SQRT
      integer i1mach, nout
C     ..
      X = T - ETA
      IF (ABS(X).LE.1.) THEN
          L = 0
          TERM = X
          S1 = TERM
   10     S = S1
          L = L + 1
          IF (L.GT.200) GO TO 20
          TERM = X*TERM/REAL(L+1)
          S1 = S + TERM
          IF (S1.NE.S) GO TO 10
          F = T*SQRT(1.+.5*THETA*T)/ (EXP(-T)*S)
          RETURN

      ELSE
          F = T*SQRT(1.+.5*THETA*T)/ (EXP(-ETA)-EXP(-T))
          RETURN

      END IF

   20 nout = i1mach(2)
      WRITE(NOUT,FMT=9000)
      RETURN

 9000 FORMAT (1X,'exp series does not converge')
      END
