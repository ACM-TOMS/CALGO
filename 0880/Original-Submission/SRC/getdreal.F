FUNCTION GETDREAL( STRING, LIST, N )
!
USE GSTEDEFINITIONS
USE DSTEDEFINITIONS
!
!.. Scalar Arguments ..
CHARACTER( * ) :: STRING
INTEGER :: N
!
!.. Derived Data Type Argument ..
TYPE( DATA_FROM_RECORD ), TARGET :: LIST
!
!.. Function Result ..
REAL( KIND=PREC ) :: GETDREAL( N )
!
!==============================================================================!
!                                                                              !
!  Purpose:                                                                    !
!  =======                                                                     !
!                                                                              !
!  GETDREAL extracts real variables from LIST.                                 !
!                                                                              !
!  Arguments:                                                                  !
!  =========                                                                   !
!                                                                              !
!  STRING  (input) CHARACTER                                                   !
!          String to be printed in case of error.                              !
!                                                                              !
!  LIST    (input) DATA_FROM_RECORD                                            !
!          List of integers.                                                   !
!                                                                              !
!  N       (input) INTEGER                                                     !
!          Number of integers in LIST.                                         !
!                                                                              !
!==============================================================================!
!  
!.. Local Scalars ..
INTEGER :: I, IOERR
!
!.. Local Derived Data Type ..
TYPE( DATA_FROM_RECORD ), POINTER :: CURRENT
!
!.. External Subroutine ..
EXTERNAL HANDLER
!
!.. Executable Statements ......................................................
!
GETDREAL = 0
CURRENT => LIST
!
DO I = 1, N
   READ ( CURRENT%FIELD, IOSTAT=IOERR, FMT=RFORMAT ) GETDREAL( I )
   IF ( IOERR /= 0 ) CALL HANDLER( 3, STRING )
   CURRENT => CURRENT%NEXT
END DO
!
END FUNCTION GETDREAL
