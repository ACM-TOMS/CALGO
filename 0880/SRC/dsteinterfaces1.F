MODULE DSTEINTERFACES1
!
!==============================================================================!
!                                                                              !
! This module defines interfaces for basic functions.                          !
!                                                                              !
!==============================================================================!
!
!-------------------------------------------------------------------------------
INTERFACE GETDREAL
   FUNCTION GETDREAL( STRING, LIST, N )
      USE GSTEDEFINITIONS
      USE DSTEDEFINITIONS
      CHARACTER( * )                    :: STRING
      INTEGER                           :: N
      REAL( KIND=PREC )                 :: GETDREAL( N )
      TYPE( DATA_FROM_RECORD ), TARGET  :: LIST
   END FUNCTION GETDREAL
END INTERFACE 
!-------------------------------------------------------------------------------
!
END MODULE DSTEINTERFACES1
