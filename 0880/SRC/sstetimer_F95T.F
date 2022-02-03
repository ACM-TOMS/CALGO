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
!.. Scalar ..
REAL :: F95_TIME
!
!.. Intrinsic Subroutine ..
INTRINSIC CPU_TIME
!
!.. Executable Statements ......................................................
!
!  Interface to Fortran 95 function.
!  On the Cray X1, TIME = SECONDR() is also an option.
!
IF ( FLAG == 0 ) THEN
   TIME = ZERO
ELSE
   CALL CPU_TIME( F95_TIME )
   TIME = REAL( F95_TIME,PREC )
END IF
!
END FUNCTION SSTETIMER
