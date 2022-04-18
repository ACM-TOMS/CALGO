MODULE GSTEINTERFACES1
!
!==============================================================================!
!                                                                              !
! This module defines interfaces for basic functions.                          !
!                                                                              !
!==============================================================================!
!
!-------------------------------------------------------------------------------
INTERFACE GETINTGR
   FUNCTION GETINTGR( STRING, LIST, N )
      USE GSTEDEFINITIONS
      CHARACTER( * )                    :: STRING
      INTEGER                           :: N
      INTEGER                           :: GETINTGR( N )
      TYPE( DATA_FROM_RECORD ), TARGET  :: LIST
   END FUNCTION GETINTGR
END INTERFACE 
!-------------------------------------------------------------------------------
INTERFACE GETSTRNG
   FUNCTION GETSTRNG( STRING, LIST, N )
      USE GSTEDEFINITIONS
      CHARACTER( LEN=FIELD_LENGTH )     :: GETSTRNG( N )
      CHARACTER( * )                    :: STRING
      INTEGER                           :: N
      TYPE( DATA_FROM_RECORD ), TARGET  :: LIST
   END FUNCTION GETSTRNG
END INTERFACE 
!-------------------------------------------------------------------------------
INTERFACE PARSER
   FUNCTION PARSER( RECORD ) RESULT( OUTPUT_FROM_PARSER )
      USE GSTEDEFINITIONS
      CHARACTER( * )                    :: RECORD
      TYPE( DATA_FROM_RECORD ), POINTER :: OUTPUT_FROM_PARSER
   END FUNCTION PARSER
END INTERFACE 
!-------------------------------------------------------------------------------
!
END MODULE GSTEINTERFACES1
