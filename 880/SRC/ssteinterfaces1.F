MODULE SSTEINTERFACES1
!
!==============================================================================!
!                                                                              !
! This module defines interfaces for basic functions.                          !
!                                                                              !
!==============================================================================!
!
!-------------------------------------------------------------------------------
INTERFACE GETSREAL
   FUNCTION GETSREAL( STRING, LIST, N )
      USE GSTEDEFINITIONS
      USE SSTEDEFINITIONS
      CHARACTER( * )                    :: STRING
      INTEGER                           :: N
      REAL( KIND=PREC )                 :: GETSREAL( N )
      TYPE( DATA_FROM_RECORD ), TARGET  :: LIST
   END FUNCTION GETSREAL
END INTERFACE 
!-------------------------------------------------------------------------------
!
END MODULE SSTEINTERFACES1
