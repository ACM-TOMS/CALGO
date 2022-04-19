c
c
      PROGRAM TEST3
c
c This program computes the n-point rational Gauss quadrature
c approximations, n=1,2,...,nmax, to the integral in Example 4.2
c for  om=.1, .5, 1., 1.999, 2., 3., 3.001, 10.5, and 50., using
c rational functions that share with the integrand the m poles
c closest to, and to the left of, the interval [0,1].
c
C     .. Parameters ..
      INTEGER NMAX,NCAPM
      PARAMETER (NMAX=10,NCAPM=50)
C     ..
C     .. Local Scalars ..
      REAL CONST,EPS,ERR,OM,SUM
      INTEGER IERAB,IERGQ,IERR,IOM,IROUT,K,KOUNT,M,MU,N,NCAP
C     ..
C     .. Local Arrays ..
      REAL A(NCAPM),ALPHA(NMAX+1),B(NCAPM),BE(NMAX+1),BETA(NMAX+1),
     +     E(NCAPM),EXACT(9),OOM(9),P0(NCAPM),P1(NCAPM),P2(NCAPM),
     +     W(NCAPM),WG(NMAX),X(NCAPM),XII(2*NMAX),XIR(2*NMAX),ZG(NMAX)
      INTEGER IS(2*NMAX)
C     ..
C     .. External Functions ..
      REAL F,R1MACH
      EXTERNAL F,R1MACH
C     ..
C     .. External Subroutines ..
      EXTERNAL ABMOD,GQRAT,RECUR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,REAL,SQRT
      integer i1mach, nout
C     ..
C     .. Data statements ..
      DATA OOM/.1,.5,1.,1.999,2.,3.,3.001,10.5,50./
      DATA EXACT/7.6724626,2.5531372,1.4784672,.81771418,.81735217,
     +     .56726445,.56709145,.17313057,.037222082/
C     ..
      nout = i1mach(2)
      WRITE(NOUT,FMT=9000)
      IROUT = 1
c      irout=0
      EPS = R1MACH(3)*1.e2
      CALL RECUR(NCAPM,6,0.,-.5,A,B,IERR)
      DO 10 K = 1,NCAPM
          A(K) = .5* (1.+A(K))
          B(K) = .25*B(K)
   10 CONTINUE
      B(1) = SQRT(8.)*B(1)
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
                          XIR(MU) = 1./OM

                      ELSE
                          XIR(MU) = 1./REAL(MU-1)
                      END IF

                      XII(MU) = 0.
                      IS(MU) = 1
   20             CONTINUE
              END IF

              CALL ABMOD(N+1,NCAPM,M,EPS,IROUT,A,B,XIR,XII,IS,ALPHA,
     +                   BETA,NCAP,KOUNT,IERAB,BE,X,W,E,P0,P1,P2)
              IF (IERAB.NE.0) WRITE(NOUT,FMT=9020) IERAB
              CALL GQRAT(N,M,ALPHA,BETA,XIR,XII,IS,ZG,WG,CONST,IERGQ,E)
              IF (IERGQ.NE.0) WRITE(NOUT,FMT=9030) IERGQ
              SUM = 0.
              DO 30 K = 1,N
                  SUM = SUM + WG(K)*F(ZG(K),OM)
   30         CONTINUE
              ERR = ABS((SUM-EXACT(IOM))/EXACT(IOM))
              IF (N.EQ.1) THEN
                  WRITE(NOUT,FMT=9040) N,M,SUM,ERR,CONST,OM

              ELSE
                  WRITE(NOUT,FMT=9050) N,M,SUM,ERR,CONST
              END IF

   40     CONTINUE
          WRITE(NOUT,FMT=9000)
   50 CONTINUE
      STOP

 9000 FORMAT (/)
 9010 FORMAT (5X,'n',4X,'m',6X,'integral',7X,'rel error',5X,'err const',
     +       /)
 9020 FORMAT (1X,'ierab=',I3)
 9030 FORMAT (1X,'iergq=',I3)
 9040 FORMAT (1X,2I5,E17.7,2E14.4,'  om=',F6.3)
 9050 FORMAT (1X,2I5,E17.7,2E14.4)
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
