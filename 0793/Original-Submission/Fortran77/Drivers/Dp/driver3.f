c
c
      PROGRAM DTEST3
c
C     .. Parameters ..
      INTEGER NMAX,NCAPM
      PARAMETER (NMAX=10,NCAPM=100)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION CONST,EPS,OM,SUM
      REAL ERR
      INTEGER IERAB,IERGQ,IERR,IOM,IROUT,K,KOUNT,M,MU,N,NCAP
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION A(NCAPM),ALPHA(NMAX+1),B(NCAPM),BE(NMAX+1),
     +                 BETA(NMAX+1),E(NCAPM),EXACT(9),OOM(9),P0(NCAPM),
     +                 P1(NCAPM),P2(NCAPM),W(NCAPM),WG(NMAX),X(NCAPM),
     +                 XII(2*NMAX),XIR(2*NMAX),ZG(NMAX)
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
      INTRINSIC ABS,DBLE,REAL,SQRT
      integer i1mach, nout
C     ..
C     .. Data statements ..
      DATA OOM/.1D0,.5D0,1.D0,1.999D0,2.D0,3.D0,3.001D0,10.5D0,50.D0/
      DATA EXACT/7.6724625863706D0,2.5531371574419D0,1.4784672063106D0,
     +     .81771417926000D0,.81735216597811D0,.56726444593607D0,
     +     .56709144671961D0,.17313056604184D0,.037222082318054D0/
C     ..
      nout = i1mach(2)
      WRITE(NOUT,FMT=9000)
      IROUT = 1
c      irout=0
      EPS = D1MACH(3)*1.D2
      CALL DRECUR(NCAPM,6,0.D0,-.5D0,A,B,IERR)
      DO 10 K = 1,NCAPM
          A(K) = .5D0* (1.D0+A(K))
          B(K) = .25D0*B(K)
   10 CONTINUE
      B(1) = SQRT(8.D0)*B(1)
      WRITE(NOUT,FMT=9010)
      DO 50 IOM = 1,9
          OM = OOM(IOM)
          DO 40 N = 1,NMAX
              M = 2*N
c          m=n
c          m=2
c          m=0
              IF (M.GT.0) THEN
                  DO 20 MU = 1,M
                      IF (MU.EQ.1) THEN
                          XIR(MU) = 1.D0/OM

                      ELSE
                          XIR(MU) = 1.D0/DBLE(MU-1)
                      END IF

                      XII(MU) = 0.D0
                      IS(MU) = 1
   20             CONTINUE
              END IF

              CALL DABMOD(N+1,NCAPM,M,EPS,IROUT,A,B,XIR,XII,IS,ALPHA,
     +                    BETA,NCAP,KOUNT,IERAB,BE,X,W,E,P0,P1,P2)
              IF (IERAB.NE.0) WRITE(NOUT,FMT=9020) IERAB
              CALL DGQRAT(N,M,ALPHA,BETA,XIR,XII,IS,ZG,WG,CONST,IERGQ,E)
              IF (IERGQ.NE.0) WRITE(NOUT,FMT=9030) IERGQ
              SUM = 0.D0
              DO 30 K = 1,N
                  SUM = SUM + WG(K)*F(ZG(K),OM)
   30         CONTINUE
              ERR = ABS(REAL((SUM-EXACT(IOM))/EXACT(IOM)))
              IF (N.EQ.1) THEN
                  WRITE(NOUT,FMT=9040) N,M,SUM,ERR,CONST,REAL(OM)

              ELSE
                  WRITE(NOUT,FMT=9050) N,M,SUM,ERR,CONST
              END IF

   40     CONTINUE
          WRITE(NOUT,FMT=9000)
   50 CONTINUE
      STOP

 9000 FORMAT (/)
 9010 FORMAT (5X,'n',4X,'m',8X,'integral',10X,'rel error',4X,
     +       'err const',/)
 9020 FORMAT (1X,'ierab=',I3)
 9030 FORMAT (1X,'iergq=',I3)
 9040 FORMAT (1X,2I5,D23.14,E13.4,D13.4,'  om=',F6.3)
 9050 FORMAT (1X,2I5,D23.14,E13.4,D13.4)
      END

      DOUBLE PRECISION FUNCTION F(T,OM)
C     .. Scalar Arguments ..
      DOUBLE PRECISION OM,T
C     ..
C     .. Local Scalars ..
      INTEGER IERR
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DGAMMA
      EXTERNAL DGAMMA
C     ..
      F = DGAMMA(1.D0+T,IERR)/ (T+OM)
      RETURN

      END
