c
c
      PROGRAM TEST5
c
c This program computes the n-point rational Gauss quadrature
c approximations, n=1,2,...,nmax, to the Fermi-Dirac integral
c using rational functions that share with the integrand the m
c conjugate complex poles closest to the interval [0,oo].
c
c      parameter(eta=1.)
c      parameter(ak=1.5)
c      parameter(ak=2.5)
C     .. Parameters ..
      INTEGER NMAX,NCAPM
      REAL THETA
      PARAMETER (NMAX=10,NCAPM=100,THETA=1.e-4)
      REAL ETA
      PARAMETER (ETA=-1.)
      REAL AK
      PARAMETER (AK=.5)
C     ..
C     .. Local Scalars ..
      REAL CONST,DEN,EPS,ERR,EXACT,PI,SUM
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
      INTRINSIC ABS,ATAN,REAL
      integer i1mach, nout
      nout = i1mach(2)
C     ..
      IROUT = 1
c      irout=0
      EPS = R1MACH(3)*1.e2
      CALL RECUR(NCAPM,7,AK,0.,A,B,IERR)
      PI = 4.*ATAN(1.)
      WRITE(NOUT,FMT=9000)
      DO 30 N = 1,NMAX
          M = 2*N
c        m=2*((n+1)/2)
c        m=2
c        m=0
          IF (M.GT.0) THEN
              DO 10 MU = 1,M - 1,2
                  DEN = ETA**2 + REAL(MU**2)* (PI**2)
                  XIR(MU) = -ETA/DEN
                  XIR(MU+1) = XIR(MU)
                  XII(MU) = REAL(MU)*PI/DEN
                  XII(MU+1) = -XII(MU)
                  IS(MU) = 1
                  IS(MU+1) = 1
   10         CONTINUE
          END IF

          CALL ABMOD(N+1,NCAPM,M,EPS,IROUT,A,B,XIR,XII,IS,ALPHA,BETA,
     +               NCAP,KOUNT,IERAB,BE,X,W,E,P0,P1,P2)
          IF (IERAB.NE.0) WRITE(NOUT,FMT=9010) IERAB
          CALL GQRAT(N,M,ALPHA,BETA,XIR,XII,IS,ZG,WG,CONST,IERGQ,E)
          IF (IERGQ.NE.0) WRITE(NOUT,FMT=9020) IERGQ
          SUM = 0.
          DO 20 K = 1,N
              SUM = SUM + WG(K)*F(ZG(K),ETA,THETA)
   20     CONTINUE
          EXACT = .29051242
c        exact=.46087845
c        exact=1.18607350
c        exact=1.39644182
c        exact=2.66187328
c        exact=7.62725610
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

      REAL FUNCTION F(T,ETA,THETA)
C     .. Scalar Arguments ..
      REAL ETA,T,THETA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC EXP,SQRT
C     ..
      F = SQRT(1.+.5*THETA*T)/ (EXP(-ETA)+EXP(-T))
      RETURN

      END
