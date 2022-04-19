FUNCTION DSTETIMER( FLAG ) RESULT( TIME )
!
USE DSTEDEFINITIONS
INCLUDE 'f90papi.h'
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
!  DSTETIMER() is the interface for a timing function.                         !
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
!.. Scalars ..
INTEGER*4 :: PAPI_INFO
INTEGER*8 :: PAPI_UTIME
!
!.. External Subroutine ..
EXTERNAL PAPIF_GET_REAL_USEC
!
!.. Executable Statements ......................................................
!
!  Interface to PAPI counters.
!
IF ( FLAG == 0 ) THEN
   TIME = ZERO
   PAPI_INFO = PAPI_VER_CURRENT
   CALL PAPIF_LIBRARY_INIT( PAPI_INFO )
   IF ( PAPI_INFO .NE. PAPI_VER_CURRENT ) THEN
      WRITE( UNIT=* , FMT='(A)' ) '* PAPI library version mismatch *'
      STOP
   END IF
ELSE
   CALL PAPIF_GET_REAL_USEC( PAPI_UTIME )
   TIME = REAL( PAPI_UTIME,PREC )*1.0D-6
END IF 
!
END FUNCTION DSTETIMER
