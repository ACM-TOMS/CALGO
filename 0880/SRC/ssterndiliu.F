SUBROUTINE SSTERNDILIU( N, ISEED, NRILIU, ILIU )
!
USE SSTEDEFINITIONS
!
!.. Scalar Arguments ..
INTEGER :: N, NRILIU
!
!.. Array Arguments ..
INTEGER :: ILIU( 2,NRILIU ), ISEED( 4 )
!
!==============================================================================!
!                                                                              !
!  Purpose:                                                                    !
!  =======                                                                     !
!                                                                              !
!  SSTERNDILIU generates random indices of the smallest and largest            !
!  eigenvalues to be computed.                                                 !
!                                                                              !
!  Arguments:                                                                  !
!  =========                                                                   !
!                                                                              !
!  N       (input) INTEGER                                                     !
!          Dimension of the matrix.                                            !
!                                                                              !
!  ISEED   (input) INTEGER array, dimension ( 4 )                              !
!          Initial seed of the random number generator. Each entry of ISEED    !
!          should lie between 0 and 4095 and ISEED(4) should be odd.           !
!                                                                              !
!  NRILIU  (input) INTEGER                                                     !
!          Number of random pair of indeces to be generated.                   !
!                                                                              !
!  ILIU    (output) INTEGER array, dimension ( 2,NRILIU )                      !
!          Indices (in ascending order) of the smallest and largest            !
!          eigenvalues to be computed, used only when RANGE='I'.               !
!                                                                              !
!==============================================================================!
!
!.. Local Scalars ..
INTEGER :: I, IL, IU
!
!.. External Function ..
REAL( KIND=PREC ), EXTERNAL :: SLARND
!
!.. Intrinsic Function ..
INTRINSIC INT
!
!.. Executable Statements ......................................................
!
DO I = 1, NRILIU
   IL = 1 + INT( ( N-1 )*SLARND( 1, ISEED ) )
   IU = 1 + INT( ( N-1 )*SLARND( 1, ISEED ) )
   IF ( IU < IL ) THEN
      ILIU( 1,I ) = IU
      ILIU( 2,I ) = IL
   ELSE
      ILIU( 1,I ) = IL
      ILIU( 2,I ) = IU
   END IF
END DO
!
END SUBROUTINE SSTERNDILIU
