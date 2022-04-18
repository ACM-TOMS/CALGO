SUBROUTINE SSTERNDVLVU( N, W, NRVLVU, VLVU, ILIU )
!
USE SSTEDEFINITIONS
!
!.. Scalar Arguments ..
INTEGER :: N, NRVLVU
!
!.. Array Arguments ..
INTEGER :: ILIU( 2,NRVLVU )
REAL( KIND=PREC ) :: W( N ), VLVU( 2,NRVLVU )
!
!==============================================================================!
!                                                                              !
!  Purpose:                                                                    !
!  =======                                                                     !
!                                                                              !
!  SSTERNDVLVU generates random bounds of intervals to be searched             !
!  for eigenvalues.                                                            !
!                                                                              !
!  Arguments:                                                                  !
!  =========                                                                   !
!                                                                              !
!  N       (input) INTEGER                                                     !
!          Dimension of the test matrix.                                       !
!                                                                              !
!  W       (input) REAL( KIND=PREC ) array, dimension ( N )                    !
!          Eigenvalues of the test matrix.                                     !
!                                                                              !
!  NRVLVU  (input) INTEGER                                                     !
!          Number of random pair of indeces to be generated.                   !
!                                                                              !
!  VLVU    (output) REAL( KIND=PREC ) array, dimension ( 2,NRVLVU )            !
!          Lower and upper bounds of intervals to be searched for              !
!          eigenvalues, used only when RANGE='V'.                              !
!                                                                              !
!  ILIU    (input) INTEGER array, dimension ( 2,NRVLVU )                       !
!          Indices (in ascending order) of the smallest and largest            !
!          eigenvalues to be computed, used only when RANGE='I'.               !
!                                                                              !
!==============================================================================!
!
!.. Local Scalars ..
INTEGER :: I, IL, IU 
REAL( KIND=PREC ) :: RTUNFL, TNORM, ULP, UNFL, VL, VU
!
!.. External Function ..
REAL( KIND=PREC ), EXTERNAL :: SLAMCH
!
!.. Intrinsic Functions ..
INTRINSIC ABS, MAX
!
!.. Executable Statements ......................................................
!
ULP = SLAMCH( 'Precision' )
UNFL = SLAMCH( 'Safe minimum' )
RTUNFL = SQRT( UNFL )
TNORM = MAX( ABS( W( 1 ) ), ABS( W( N ) ) )
!
DO I = 1, NRVLVU
   IL = ILIU( 1,I )
   IU = ILIU( 2,I )
   IF ( IL /= 1 ) THEN
      VL = W( IL ) - MAX( HALF*( W( IL )-W( IL-1 ) ), &
                          ULP*TNORM, TWO*RTUNFL )
   ELSE
      VL = W( 1 ) - MAX( HALF*( W( N )-W( 1 ) ), &
                         ULP*TNORM, TWO*RTUNFL )
   END IF
   IF ( IU /= N ) THEN
      VU = W( IU ) + MAX( HALF*( W( IU+1 )-W( IU ) ), &
                          ULP*TNORM, TWO*RTUNFL )
   ELSE
      VU = W( N ) + MAX( HALF*( W( N )-W( 1 ) ), &
                         ULP*TNORM, TWO*RTUNFL )
   END IF
   VLVU( 1,I ) = VL
   VLVU( 2,I ) = VU
END DO
!
END SUBROUTINE SSTERNDVLVU
