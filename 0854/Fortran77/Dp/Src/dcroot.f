      SUBROUTINE DCROOT( XR, XI, YR, YI )
      IMPLICIT NONE
C
C     PURPOSE
C
C     To compute the complex square root YR + i*YI of a complex number
C     XR + i*XI  in real arithmetic.  The branch is chosen so that
C     YR .GE. 0.0  and  SIGN(YI) .EQ. SIGN(XI).
C
C     ARGUMENTS
C
C     XR      (input) DOUBLE PRECISION
C     XI      (input) DOUBLE PRECISION
C             On entry, these scalars define the real and imaginary
C             part of the complex number of which the root is to be
C             computed.
C
C     YR      (output) DOUBLE PRECISION
C     YI      (output) DOUBLE PRECISION
C             These scalars define the real and imaginary part of the
C             complex root.
C
C     METHOD
C
C     This routine is identical (apart from the documentation) with
C     the routine DCROOT in Algorithm 800, TOMS.
C
C     CONTRIBUTORS
C
C     D. Kressner, Technical Univ. Berlin, Germany, and
C     P. Benner, Technical Univ. Chemnitz, Germany, December 2003.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO,HALF
      PARAMETER          ( ZERO = 0.0D0, HALF = 0.5D0 )
C     .. Scalar Arguments ..
      DOUBLE PRECISION   XI, XR, YI, YR
C     .. Local Scalars ..
      DOUBLE PRECISION   S
C     .. External Functions ..
      DOUBLE PRECISION   DLAPY2
      EXTERNAL           DLAPY2

C     .. Intrinsic Functions ..
      INTRINSIC          ABS, SQRT
C
C     .. Executable Statements ..
C
      S = SQRT( HALF*( DLAPY2( XR, XI ) + ABS( XR ) ) )
      IF ( XR.GE.ZERO )
     $   YR = S
      IF ( XI.LT.ZERO )
     $   S = -S
      IF ( XR.LE.ZERO )
     $   YI = S
      IF ( XR.LT.ZERO )
     $   YR = HALF* (XI/YI)
      IF ( XR.GT.ZERO )
     $   YI = HALF* (XI/YR)
C
      RETURN
C *** Last line of DCROOT ***
      END
