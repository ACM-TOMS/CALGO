      PROGRAM DRIVER
C     .. Parameters ..
      INTEGER NMAX
      PARAMETER (NMAX=100)
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION A,B,IN,IVAL
      INTEGER I,N,P,R
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION C(0:NMAX+1),T(0:NMAX),W(0:NMAX),WSP(0:NMAX)
C     ..
C     .. External Subroutines ..
      EXTERNAL GREGORY
C     ..
      A = -ONE
      B = ONE
      N = 7
      R = N
C
C limits of integration
C
      T(0) = A
      T(N) = B
C
C Generate the weights
C
      CALL GREGORY(N,R,T,W,WSP,C)
      DO I = 0,N
          WRITE (*,FMT='(2e14.6)') W(I),T(I)
      END DO
C
C Check - integrate x^p */
C
      DO P = 0,R + 4
          IVAL = (B** (P+1)-A** (P+1))/ (P+1)
          IN = 0
          DO I = 0,N
              IN = IN + W(I)*T(I)**P
          END DO
          WRITE (*,FMT='(3i4, 2e14.6, e12.4)') N,R,P,IVAL,IN,IVAL - IN
      END DO
      END
