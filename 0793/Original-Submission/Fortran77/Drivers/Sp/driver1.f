c
c
      PROGRAM TEST1
c
c This program computes the n-point rational Gauss quadrature
c approximations, n=1,2,...,nmax, to the integral in Example 4.1
c using rational functions that share with the integrand the m
c poles closest to the interval [-1,1].
c
C     .. Parameters ..
      INTEGER NMAX,NCAPM
      REAL OM
      PARAMETER (NMAX=10,NCAPM=50,OM=1.1)
C     ..
C     .. Local Scalars ..
      REAL CONST,EPS,ERR,EXACT,SGN,SUM
      INTEGER IERAB,IERGQ,IERR,IROUT,K,KOUNT,M,MU,N,NCAP
C     ..
C     .. Local Arrays ..
      REAL A(NCAPM),ALPHA(NMAX+1),B(NCAPM),BE(NMAX+1),BETA(NMAX+1),
     +     E(NCAPM),P0(NCAPM),P1(NCAPM),P2(NCAPM),W(NCAPM),WG(NMAX),
     +     X(NCAPM),XII(2*NMAX),XIR(2*NMAX),ZG(NMAX)
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
      INTRINSIC ABS,REAL
      integer i1mach, nout
      nout = i1mach(2)
C     ..
      IROUT = 1
c      irout=0
      EPS = R1MACH(3)*1.e2
      CALL RECUR(NCAPM,1,0.,0.,A,B,IERR)
      WRITE(NOUT,FMT=9000)
      DO 30 N = 1,NMAX
          M = 2*N
c        m=2*((n+1)/2)
c        m=2
c        m=0
          SGN = 1.
          DO 10 MU = 1,M
              SGN = -SGN
              XIR(MU) = SGN/ (OM*REAL((MU+1)/2))
              XII(MU) = 0.
              IS(MU) = 1
   10     CONTINUE
          CALL ABMOD(N+1,NCAPM,M,EPS,IROUT,A,B,XIR,XII,IS,ALPHA,BETA,
     +               NCAP,KOUNT,IERAB,BE,X,W,E,P0,P1,P2)
          IF (IERAB.NE.0) WRITE(NOUT,FMT=9010) IERAB
          CALL GQRAT(N,M,ALPHA,BETA,XIR,XII,IS,ZG,WG,CONST,IERGQ,E)
          IF (IERGQ.NE.0) WRITE(NOUT,FMT=9020) IERGQ
          SUM = 0.
          DO 20 K = 1,N
              SUM = SUM + WG(K)*F(ZG(K),OM)
   20     CONTINUE
          EXACT = 4.4677736
          ERR = ABS((SUM-EXACT)/EXACT)
          WRITE(NOUT,FMT=9030) N,M,SUM,ERR,CONST
   30 CONTINUE
      STOP

 9000 FORMAT (5X,'n',4X,'m',6X,'integral',7X,'rel error',5X,'err const',
     +       /)
 9010 FORMAT (1X,'ierab=',I3)
 9020 FORMAT (1X,'iergq=',I3)
 9030 FORMAT (1X,2I5,E17.7,2E14.4)
      END

      REAL FUNCTION F(T,OM)
C     .. Scalar Arguments ..
      REAL OM,T
C     ..
C     .. Local Scalars ..
      REAL PI
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ATAN,SIN
C     ..
      PI = 4.*ATAN(1.)
      F = 1.
      IF (T.NE.0.) F = (PI*T/OM)/SIN(PI*T/OM)
      RETURN

      END
