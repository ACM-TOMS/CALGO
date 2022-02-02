c
c
      PROGRAM DTEST1
c
c This is a double-precision version of test1
c
C     .. Parameters ..
      INTEGER NMAX,NCAPM
      PARAMETER (NMAX=10,NCAPM=100)
      DOUBLE PRECISION OM
      PARAMETER (OM=1.1D0)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION CONST,EPS,EXACT,SGN,SUM
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
      INTRINSIC ABS,DBLE,REAL
      integer i1mach, nout
      nout = i1mach(2)
C     ..
      IROUT = 1
c      irout=0
      EPS = D1MACH(3)*1.D2
      CALL DRECUR(NCAPM,1,0.D0,0.D0,A,B,IERR)
      WRITE(NOUT,FMT=9000)
      DO 30 N = 1,NMAX
          M = 2*N
c        m=n
c        m=2
c        m=0
          SGN = 1.D0
          DO 10 MU = 1,M
              SGN = -SGN
              XIR(MU) = SGN/ (OM*DBLE((MU+1)/2))
              XII(MU) = 0.D0
              IS(MU) = 1
   10     CONTINUE
          CALL DABMOD(N+1,NCAPM,M,EPS,IROUT,A,B,XIR,XII,IS,ALPHA,BETA,
     +                NCAP,KOUNT,IERAB,BE,X,W,E,P0,P1,P2)
          IF (IERAB.NE.0) WRITE(NOUT,FMT=9010) IERAB
          CALL DGQRAT(N,M,ALPHA,BETA,XIR,XII,IS,ZG,WG,CONST,IERGQ,E)
          IF (IERGQ.NE.0) WRITE(NOUT,FMT=9020) IERGQ
          SUM = 0.D0
          DO 20 K = 1,N
              SUM = SUM + WG(K)*F(ZG(K),OM)
   20     CONTINUE
          EXACT = 4.467773646387766D0
          ERR = ABS(REAL((SUM-EXACT)/EXACT))
          WRITE(NOUT,FMT=9030) N,M,SUM,ERR,CONST
   30 CONTINUE
      STOP

 9000 FORMAT (5X,'n',4X,'m',8X,'integral',9X,'rel error',3X,'err const',
     +       /)
 9010 FORMAT (1X,'ierab=',I3)
 9020 FORMAT (1X,'iergq=',I3)
 9030 FORMAT (1X,2I5,D23.14,E12.4,D12.4)
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
