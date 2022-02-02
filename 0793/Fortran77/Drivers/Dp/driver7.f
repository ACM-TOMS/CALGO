c
c
      PROGRAM DTEST7
c
c This is a double-precision version of  test7
c
c
c Be sure that ncapm is greater or equalt to numax
c
c      parameter(ak=1.5d0)
c      parameter(ak=2.5d0)
C     .. Parameters ..
      INTEGER NMAX,NCAPM,NUMAX
      PARAMETER (NMAX=10,NCAPM=100,NUMAX=50)
      DOUBLE PRECISION THETA
      PARAMETER (THETA=1.D-4)
      DOUBLE PRECISION AK
      PARAMETER (AK=.5D0)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION C1,C2,CONST,DX,DY,EPS,ETA,EXACT,FMU,FMU2,FNU,
     +                 FNU2,HN,HP,P,PI,SQPI,SUM,SUM0,T,X1,XA
      REAL ERR
      INTEGER IERAB,IERGC,IERGQ,IERK,IERR,IROUT,J,K,KOUNT,M,MU,N,NCAP,
     +        NU,NU0,NU1
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION A(NCAPM),A1(NMAX+1),ALPHA(NMAX+1),B(NCAPM),
     +                 B1(NMAX+1),BE(NMAX+1),BETA(NMAX+1),E(NCAPM),
     +                 P0(NCAPM),P1(NCAPM),P2(NCAPM),RHOI(1),RHOR(1),
     +                 ROLDI(1),ROLDR(1),W(NCAPM),WG(NMAX),X(NCAPM),
     +                 XII(2*NMAX),XIR(2*NMAX),ZG(NMAX)
      INTEGER IS(2*NMAX)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION D1MACH,F
      EXTERNAL D1MACH,F
C     ..
C     .. External Subroutines ..
      EXTERNAL DABMOD,DGCHRS,DGQRAT,DKNUM,DRECUR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,ATAN,DBLE,EXP,REAL,SQRT
      integer i1mach, nout
      nout = i1mach(2)
C     ..
      PI = 4.D0*ATAN(1.D0)
      SQPI = SQRT(PI)
      C1 = PI
c      c1=-pi
c      c1=pi
      C2 = SQPI
c      c2=.5*sqpi
c      c2=.75*sqpi
      IROUT = 1
c      irout=0
      EPS = D1MACH(3)*1.D3
      CALL DRECUR(NCAPM,7,AK-1.D0,0.D0,A,B,IERR)
      WRITE(NOUT,FMT=9000)
      DO 60 N = 1,NMAX
          M = 2*N - 1
c        m=2*((n+1)/2)-1
c        m=1
          IF (M.GE.3) THEN
              DO 10 MU = 1,M - 2,2
                  XIR(MU) = 0.D0
                  XIR(MU+1) = 0.D0
                  XII(MU) = 1.D0/ (DBLE(MU+1)*PI)
                  XII(MU+1) = -XII(MU)
                  IS(MU) = 1
                  IS(MU+1) = 1
   10         CONTINUE
          END IF

          CALL DABMOD(N+1,NCAPM,M-1,EPS,IROUT,A,B,XIR,XII,IS,A1,B1,NCAP,
     +                KOUNT,IERAB,BE,X,W,E,P0,P1,P2)
          IF (IERAB.NE.0) WRITE(NOUT,FMT=9010) IERAB
          ETA = -.001D0
          XIR(M) = -1.D0/ETA
          XII(M) = 0.D0
          IS(M) = 1
          X1 = ETA
          XA = ABS(ETA)
          T = 1.D0
          SUM0 = 0.D0
          SUM = -1.D0/ (AK-1.D0)
          J = 0
   20     J = J + 1
          IF (J.GT.100) THEN
              WRITE(NOUT,FMT=9020)
              GO TO 60

          END IF

          T = -XA*T/DBLE(J)
          SUM = SUM + T/ (DBLE(J+1)-AK)
          IF (SUM.NE.SUM0) THEN
              SUM0 = SUM
              GO TO 20

          END IF

          HP = -EXP(XA)* (C1* (XA** (AK-1.D0))-C2*SUM)
          HN = 0.D0
          IF (M.GE.3) THEN
              SUM = 0.D0
              DO 40 NU = 1, (M-1)/2
                  FNU = DBLE(NU)
                  FNU2 = FNU**2
                  P = 1.D0
                  DO 30 MU = 1, (M-1)/2
                      IF (MU.EQ.NU) GO TO 30
                      FMU = DBLE(MU)
                      FMU2 = FMU**2
                      P = (FMU2/ (FMU2-FNU2))*P
   30             CONTINUE
                  NU0 = 10
                  DX = 0.D0
                  DY = 2.D0*FNU*PI
                  CALL DKNUM(0,NU0,NUMAX,DX,DY,EPS,A,B,RHOR,RHOI,NU1,
     +                       IERK,ROLDR,ROLDI)
                  IF (IERK.NE.0) WRITE(NOUT,FMT=9030) IERK
                  SUM = SUM + 2.D0*FNU*PI*P*
     +                  (2.D0*FNU*PI*HP- (2.D0*FNU*PI*RHOR(1)+
     +                  X1*RHOI(1)))/ (X1**2+4.D0*FNU2* (PI**2))
   40         CONTINUE
              HP = SUM
          END IF

          CALL DGCHRS(N+1,1,A1,B1,X1,HP,HN,ALPHA,BETA,IERGC)
          IF (IERGC.NE.0) WRITE(NOUT,FMT=9040) IERGC
          BETA(1) = -ETA*BETA(1)
          CALL DGQRAT(N,M,ALPHA,BETA,XIR,XII,IS,ZG,WG,CONST,IERGQ,E)
          IF (IERGQ.NE.0) WRITE(NOUT,FMT=9050) IERGQ
          SUM = 0.D0
          DO 50 K = 1,N
              SUM = SUM + WG(K)*F(ZG(K),ETA,THETA)
   50     CONTINUE
          EXACT = 2.2171501009112D0
c        exact=1.7800123286291d0
c        exact=3.7403844583705d0
          ERR = ABS(REAL((SUM-EXACT)/EXACT))
          WRITE(NOUT,FMT=9060) N,M,SUM,ERR,CONST
   60 CONTINUE
      STOP

 9000 FORMAT (5X,'n',4X,'m',8X,'integral',10X,'rel error',4X,
     +       'err const',/)
 9010 FORMAT (1X,'ierab=',I3)
 9020 FORMAT (1X,'power series for gamma does not converge')
 9030 FORMAT (1X,'ierk for complex z =',I3)
 9040 FORMAT (1X,'iergc=',I3)
 9050 FORMAT (1X,'iergq=',I3)
 9060 FORMAT (1X,2I5,D23.14,E13.4,D13.4)
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
