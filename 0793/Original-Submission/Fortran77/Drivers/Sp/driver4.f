c
c
      PROGRAM TEST4
c
c This program does the same as  test3  but treats the pole closest
c to the interval [0,1] separately.
c
C     .. Parameters ..
      INTEGER NMAX,NCAPM
      PARAMETER (NMAX=10,NCAPM=50)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION DHP,DOM,DP,DSUM,DX1
      REAL CONST,EPS,ERR,EXACT,HN,HP,OM,SUM,X1
      INTEGER IERAB,IERGC,IERGQ,IERR,IROUT,K,KOUNT,M,MU,N,NCAP,NU
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION DXIR(2*NMAX)
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
      EXTERNAL ABMOD,GCHRIS,GQRAT,RECUR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,ATAN,DBLE,REAL,SQRT
      integer i1mach, nout
      nout = i1mach(2)
C     ..
      IROUT = 1
c      irout=0
      EPS = R1MACH(3)*1.e2
      CALL RECUR(NCAPM,6,0.,-.5,A,B,IERR)
      DO 10 K = 1,NCAPM
          A(K) = .5* (1.+A(K))
          B(K) = .25*B(K)
   10 CONTINUE
      B(1) = SQRT(8.)*B(1)
      WRITE(NOUT,FMT=9000)
      DO 60 N = 1,NMAX
          M = 2*N
c        m=n
c        m=1
          IF (M.GT.0) THEN
              DO 20 NU = 1,M
                  DXIR(NU) = 1.D0/DBLE(NU)
                  XIR(NU) = REAL(DXIR(NU))
                  XII(NU) = 0.
                  IS(NU) = 1
   20         CONTINUE
          END IF

          CALL ABMOD(N+1,NCAPM,M-1,EPS,IROUT,A,B,XIR,XII,IS,A1,B1,NCAP,
     +               KOUNT,IERAB,BE,X,W,E,P0,P1,P2)
          IF (IERAB.NE.0) WRITE(NOUT,FMT=9010) IERAB,N
          DOM = .001D0
          OM = REAL(DOM)
          DXIR(M) = 1.D0/DOM
          XIR(M) = REAL(DXIR(M))
          DX1 = -DOM
          X1 = REAL(DX1)
          DHP = -2.D0*ATAN(1.D0/SQRT(ABS(DX1)))/SQRT(ABS(DX1))
          HP = REAL(DHP)
          HN = 0.
          IF (M.GT.1) THEN
              DSUM = 0.D0
              DO 40 NU = 1,M - 1
                  DP = 1.D0
                  DO 30 MU = 1,M - 1
                      IF (MU.EQ.NU) GO TO 30
                      DP = (1.D0-DXIR(MU)/DXIR(NU))*DP
   30             CONTINUE
                  DSUM = DSUM + (DHP+2.D0*SQRT(DXIR(NU))*
     +                   ATAN(SQRT(DXIR(NU))))/ ((1.D0+DXIR(NU)*DX1)*DP)
   40         CONTINUE
              DHP = DSUM
              HP = REAL(DHP)
          END IF

          CALL GCHRIS(N+1,1,A1,B1,X1,HP,HN,ALPHA,BETA,IERGC)
          IF (IERGC.NE.0) WRITE(NOUT,FMT=9020) IERGC,N
          BETA(1) = OM*BETA(1)
          CALL GQRAT(N,M,ALPHA,BETA,XIR,XII,IS,ZG,WG,CONST,IERGQ,E)
          IF (IERGQ.NE.0) WRITE(NOUT,FMT=9030) IERGQ,N
          SUM = 0.
          DO 50 K = 1,N
              SUM = SUM + WG(K)*F(ZG(K),OM)
   50     CONTINUE
          EXACT = 96.703688
          ERR = ABS((SUM-EXACT)/EXACT)
          WRITE(NOUT,FMT=9040) N,M,SUM,ERR,CONST
   60 CONTINUE
      STOP

 9000 FORMAT (5X,'n',4X,'m',6X,'integral',7X,'rel error',5X,'err const',
     +       /)
 9010 FORMAT (1X,'ierab=',I3,'  for n=',I2)
 9020 FORMAT (1X,'iergc=',I3,'  for n=',I2)
 9030 FORMAT (1X,'iergq=',I3,'  for n=',I2)
 9040 FORMAT (1X,2I5,E17.7,2E14.4)
      END

      REAL FUNCTION F(T,OM)
C     .. Scalar Arguments ..
      REAL OM,T
C     ..
C     .. Local Scalars ..
      INTEGER IERR
C     ..
C     .. External Functions ..
      REAL GAMMA
      EXTERNAL GAMMA
C     ..
      F = GAMMA(1.+T,IERR)/ (T+OM)
      RETURN

      END
