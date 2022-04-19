      PROGRAM MAIN
C Driver for  CACM Alg 125 Rutishauser
C Author M.Dow@anu.edu.au
C        ANUSF,  Australian National University
Canberra Australia
C
C Tidied up to use workspace arrays, real parameter values
C and put through nag tools
C
C trh (20/07/97)
C
C     ..
C     .. Parameters ..
      INTEGER NM
      PARAMETER (NM=100)
      DOUBLE PRECISION ZERO,ONE,TWO,THREE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,THREE=3.0D0)
C     .. Local Scalars ..
      DOUBLE PRECISION A,B,EPS,EXACT,S
      INTEGER I,N,P
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION E(NM),Q(NM),W(NM),WORK(9*NM+8),X(NM)
C     ..
C     .. External Subroutines ..
      EXTERNAL WEIGHTCOEFF
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MIN
C     ..
   10 READ (5,FMT=*,END=20) N,EPS
      IF (EPS.LE.ZERO) EPS = 1D-15
      IF (N.GE.NM-1) STOP 99
      WRITE (6,FMT=*) 'n=',N,' eps=',EPS
C
C  q,e for w=1/2 interval (0,2)
C Ref: W.Jones & W.Thron Continued Fractions ...
C      Encyclopedia of maths...vol 11 p24
C      Continued Fraction expansion of log(1+x),
C transformed to Rutishauser form p33
C
      DO I = 2,NM
          Q(I) = TWO*I*I/ (TWO*I* (TWO*I-ONE))
          E(I) = TWO*I*I/ (TWO*I* (TWO*I+ONE))
      END DO
      Q(1) = ONE
      E(1) = ONE/THREE
      CALL WEIGHTCOEFF(N,Q,E,EPS,W,X,WORK)
C Adjust weights, zeros for w=1, interval (-1,1)
      DO I = 1,N
          W(I) = TWO*W(I)
          X(I) = X(I) - ONE
      END DO

      DO I = 1,MIN(N,10)
          WRITE (6,FMT=9000) W(I),X(I)
      END DO
C Check for x^4
      P = 4
      A = -ONE
      B = ONE
      EXACT = (B** (P+1)-A** (P+1))/ (P+1)
      S = ZERO
      DO I = 1,N
          S = S + W(I)*X(I)**P
      END DO
      WRITE (6,FMT=9010) P,EXACT,S,EXACT - S
      GO TO 10

   20 STOP

 9000 FORMAT (F20.16,2X,F20.16)
 9010 FORMAT (I3,' Exact=',F20.16,' Quadrature=',F20.16,' Error=',E14.7)
      END
