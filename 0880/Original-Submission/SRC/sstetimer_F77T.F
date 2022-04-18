FUNCTION SSTETIMER( FLAG ) RESULT( TIME )
!
USE SSTEDEFINITIONS
!
!.. Scalar Argument ..
INTEGER :: FLAG  
!
!.. Function Result ..
REAL( KIND=PREC ) :: TIME
!
!===============================================================================
!                                                                              !
!  Purpose:                                                                    !
!  =======                                                                     !
!                                                                              !
!  SSTETIMER() is the interface for a timing function.                         !
!                                                                              !
!  Argument:                                                                   !
!  ========                                                                    !
!                                                                              !
!  FLAG  (input) INTEGER                                                       !
!        = 0, initialize timer (required for PAPI)                             !
!        = 1, get time                                                         !
!                                                                              !
!===============================================================================
!  
!.. Static Array ..
REAL :: TARRAY( 2 )
!
!.. External Function ..
EXTERNAL ETIME  
!
!.. Executable Statements ......................................................
!
!  Interface to standard Fortran functions.
!
IF ( FLAG == 0 ) THEN
   TIME = ZERO
ELSE
   TIME = ETIME( TARRAY )
   TIME = REAL( TARRAY( 1 ),PREC )
END IF
!
END FUNCTION SSTETIMER
