c
c
      PROGRAM DTEST2
c
c This is a double-precision version of  test2
c
C     .. Parameters ..
      INTEGER NMAX,NCAPM
      PARAMETER (NMAX=10,NCAPM=100)
      DOUBLE PRECISION OM
      PARAMETER (OM=1.001D0)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION CONST,EPS,EXACT,HN,HP,P,SGN,SUM,SUMN,SUMP,X1
      REAL ERR
      INTEGER IERAB,IERGC,IERGQ,IERR,IROUT,K,KOUNT,M,MU,N,NCAP,NU
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION A(NCAPM),A1(NMAX+1),ALPHA(NMAX+1),B(NCAPM),
     +                 B1(NMAX+1),BE(NMAX+1),BETA(NMAX+1),E(NCAPM),
     +                 P0(NCAPM),P1(NCAPM),P2(NCAPM),W(NCAPM),WG(NMAX),
     +                 X(NCAPM),XII(2*NMAX),XIR(2*NMAX),ZG(NMAX)
      INTEGER IS(2*NMAX)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION D1MACH,F
      EXTERNAL D1MACH,F
C     ..
C     .. External Subroutines ..
      EXTERNAL DABMOD,DGCHRS,DGQRAT,DRECUR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,DBLE,LOG,REAL
      integer i1mach, nout
      nout = i1mach(2)
C     ..
      IROUT = 1
c      irout=0
      EPS = D1MACH(3)*1.D2
      CALL DRECUR(NCAPM,1,0.D0,0.D0,A,B,IERR)
      WRITE(NOUT,FMT=9000)
      DO 50 N = 1,NMAX
          M = 2*N
c        m=n
c        m=4
c        m=2
          SGN = 1.D0
          DO 10 NU = 1,M
              SGN = -SGN
              IF (NU.LE.M-2) XIR(NU) = SGN/DBLE((NU+3)/2)
              XII(NU) = 0.D0
              IS(NU) = 1
   10     CONTINUE
          CALL DABMOD(N+1,NCAPM,M-2,EPS,IROUT,A,B,XIR,XII,IS,A1,B1,NCAP,
     +                KOUNT,IERAB,BE,X,W,E,P0,P1,P2)
          IF (IERAB.NE.0) WRITE(NOUT,FMT=9010) IERAB
          XIR(M-1) = 1.D0/OM
          XIR(M) = -1.D0/OM
          X1 = OM
          HP = LOG(ABS((X1+1.D0)/ (X1-1.D0)))
          HN = -HP
          IF (M.GT.2) THEN
              SUMP = 0.D0
              SUMN = 0.D0
              DO 30 NU = 1,M - 2
                  P = 1.D0
                  DO 20 MU = 1,M - 2
                      IF (MU.EQ.NU) GO TO 20
                      P = (1.D0-XIR(MU)/XIR(NU))*P
   20             CONTINUE
                  SUMP = SUMP + LOG(ABS((X1+1.D0)* (XIR(NU)+1.D0)/ ((X1-
     +                   1.D0)* (XIR(NU)-1.D0))))/ ((1.D0+XIR(NU)*X1)*P)
                  SUMN = SUMN + LOG(ABS((X1-1.D0)* (XIR(NU)+1.D0)/ ((X1+
     +                   1.D0)* (XIR(NU)-1.D0))))/ ((1.D0-XIR(NU)*X1)*P)
   30         CONTINUE
              HP = SUMP
              HN = SUMN
          END IF

          CALL DGCHRS(N+1,2,A1,B1,X1,HP,HN,ALPHA,BETA,IERGC)
          BETA(1) = - (OM**2)*BETA(1)
          CALL DGQRAT(N,M,ALPHA,BETA,XIR,XII,IS,ZG,WG,CONST,IERGQ,E)
          IF (IERGQ.NE.0) WRITE(NOUT,FMT=9020) IERGQ
          SUM = 0.D0
          DO 40 K = 1,N
              SUM = SUM + WG(K)*F(ZG(K),OM)
   40     CONTINUE
          EXACT = 12.929256850003D0
          ERR = ABS(REAL((SUM-EXACT)/EXACT))
          WRITE(NOUT,FMT=9030) N,M,SUM,ERR,CONST
   50 CONTINUE
      STOP

 9000 FORMAT (5X,'n',4X,'m',8X,'integral',10X,'rel error',4X,
     +       'err const',/)
 9010 FORMAT (1X,'ierab=',I3)
 9020 FORMAT (1X,'iergq=',I3)
 9030 FORMAT (1X,2I5,D23.14,E13.4,D13.4)
      END

      DOUBLE PRECISION FUNCTION F(T,OM)
C     .. Scalar Arguments ..
      DOUBLE PRECISION OM,T
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION PI
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ATAN,SIN
C     ..
      PI = 4.D0*ATAN(1.D0)
      F = 1.D0
      IF (T.NE.0.D0) F = (PI*T/OM)/SIN(PI*T/OM)
      RETURN

      END
