      SUBROUTINE ZLACPV( N, I, J, A, B )
C
C     PURPOSE
C
C     To copy N entries starting at the index I of the real vector A to
C     the complex vector B starting at the index J.
C
C     This is used for calls involving %VAL constructs.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     N      (input) INTEGER
C            The number of entries to copy.  N >= 0.
C
C     I      (input) INTEGER
C            The starting index in the vector A.
C
C     J      (input) INTEGER
C            The starting index in the vector B.
C
C     A      (input) DOUBLE PRECISION array, dimension (*)
C            The given vector A.
C
C     B      (output) COMPLEX*16 array, dimension (*)
C            The resulted vector B:
C            B(J) = A(I), ..., B(J+N-1) = A(I+N-1).
C
C     CONTRIBUTOR
C
C     V. Sima, Research Institute for Informatics, Bucharest, Aug. 2012.
C
C     REVISIONS
C
C     V. Sima, Nov. 2012, June 2014.
C     M. Voigt, Jul. 2013.
C
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER            DKIND
      PARAMETER          ( DKIND = KIND( 0.0D+0 ) )
C     .. Scalar Arguments ..
      INTEGER            I, J, N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   A( * )
      COMPLEX(DKIND) ::  B( * )
C
C     .. Local Scalars ..
      INTEGER            K, L
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          CMPLX
C
C     .. Executable Statements ..
C
C     For efficiency, the input parameters are not checked.
C
      L = J
      DO 10 K = I, I+N-1
         B(L) =  CMPLX( A(K), KIND = DKIND )
         L = L + 1
   10 CONTINUE
C
      RETURN
C *** Last line of ZLACPV ***
      END
