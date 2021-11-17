MODULE GSTEINTERFACES2
!
!==============================================================================!
!                                                                              !
! This module defines interfaces for basic functions.                          !
!                                                                              !
!==============================================================================!
!
!-------------------------------------------------------------------------------
INTERFACE PARSERLIST
   FUNCTION PARSERLIST ( STRING, LIST ) RESULT( OUTPUT_FROM_PARSERLIST )
      USE GSTEDEFINITIONS
      USE GSTEINTERFACES1, ONLY : GETINTGR, PARSER
      CHARACTER( * )                    :: STRING
      TYPE( DATA_FROM_RECORD ), TARGET  :: LIST
      TYPE( DATA_FROM_RECORD ), POINTER :: OUTPUT_FROM_PARSERLIST
   END FUNCTION PARSERLIST
END INTERFACE
!-------------------------------------------------------------------------------
!
END MODULE GSTEINTERFACES2
