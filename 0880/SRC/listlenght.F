INTEGER FUNCTION LISTLENGHT( LIST ) 
!
USE GSTEDEFINITIONS
!
!.. Derived Data Type Argument ..
TYPE( DATA_FROM_RECORD ) :: LIST
!
!==============================================================================!
!                                                                              !
!  Purpose:                                                                    !
!  =======                                                                     !
!                                                                              !
!  LISTLENGHT finds out the number of items in LIST.                           !
!                                                                              !
!  Arguments:                                                                  !
!  =========                                                                   !
!                                                                              !
!  LIST  (input) DATA_FROM_RECORD                                              !
!        List of data as interpreted by PARSER.                                !
!                                                                              !
!===============================================================================
!
!.. Local Scalar ..
INTEGER :: NDATA
!
!.. Local Derived Data Type ..
TYPE( DATA_FROM_RECORD ), POINTER :: LIST_CURRENT
!
!.. Executable Statements ......................................................
!
NDATA = 1
LIST_CURRENT => LIST%NEXT
DO WHILE ( ASSOCIATED( LIST_CURRENT ) )
   NDATA = NDATA + 1
   LIST_CURRENT => LIST_CURRENT%NEXT
END DO
LISTLENGHT = NDATA
!
END FUNCTION LISTLENGHT
