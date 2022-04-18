FUNCTION GETINTGR( STRING, LIST, N )
!
USE GSTEDEFINITIONS
!
!.. Scalar Arguments ..
CHARACTER( * ) :: STRING
INTEGER :: N
!
!.. Derived Data Type Argument ..
TYPE( DATA_FROM_RECORD ), TARGET :: LIST
!
!.. Function Result ..
INTEGER :: GETINTGR( N )
!
!==============================================================================!
!                                                                              !
!  Purpose:                                                                    !
!  =======                                                                     !
!                                                                              !
!  GETINTGR extracts integer values from LIST.                                 !
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
TYPE( DATA_FROM_RECORD ), POINTER :: LIST_CURRENT
!
!.. External Subroutine ..
EXTERNAL HANDLER
!
!.. Executable Statements ......................................................
!
GETINTGR = 0
LIST_CURRENT => LIST
!
DO I = 1, N
   READ ( LIST_CURRENT%FIELD, IOSTAT=IOERR, FMT=IFORMAT ) GETINTGR( I )
   IF ( IOERR /= 0 ) CALL HANDLER( 3, STRING )
   LIST_CURRENT => LIST_CURRENT%NEXT
END DO
!
END FUNCTION GETINTGR
