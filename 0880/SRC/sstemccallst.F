SUBROUTINE SSTEMCCALLST( NTEST, MACRO, LIST, CALLST )
!
USE GSTEDEFINITIONS
USE GSTEINTERFACES1, ONLY : GETSTRNG
! 
!.. Scalar Arguments ..
CHARACTER( LEN=MACRO_NAME_LENGTH ) :: MACRO
INTEGER :: NTEST
!
!.. Array Arguments ..
LOGICAL :: CALLST( 8 )
!
!.. Derived Data Type Argument ..
TYPE( DATA_FROM_RECORD ), POINTER :: LIST
!
!==============================================================================!
!                                                                              !
!  Purpose:                                                                    !
!  =======                                                                     !
!                                                                              !
!  SSTEMCCALLST deals with the macro that sets tests to be performed.          !
!                                                                              !
!  Arguments:                                                                  !
!  =========                                                                   !
!                                                                              !
!  NTEST   (output) INTEGER                                                    !
!          Number of tests to be performed.                                    !
!                                                                              !
!  MACRO   (input) CHARACTER                                                   !
!          Macro definition.                                                   !
!                                                                              !
!  LIST    (input) DATA_FROM_RECORD (derived data type)                        !
!          List of strings.                                                    !
!                                                                              !
!  CALLST  (output) LOGICAL, dimension ( 8 )                                   !
!          Sets tests to be performed.                                         !
!                                                                              !
!==============================================================================!
! 
!.. Local Scalar ..
INTEGER :: I
!
!.. Static Array ..
CHARACTER( LEN=MACRO_NAME_LENGTH ) :: TEST_LIST( 8 )
!
!.. External Suboroutine ..
EXTERNAL HANDLER
!
!.. External Function ..
INTEGER, EXTERNAL :: LISTLENGHT
!
!.. Executable Statements ......................................................
!
NTEST = LISTLENGHT( LIST )
TEST_LIST( 1:NTEST ) = GETSTRNG( MACRO, LIST, NTEST )
DO I = 1, NTEST
   SELECT CASE ( TEST_LIST( I ) )
   CASE ( 'ALL', 'all' ); CALLST = .TRUE.; EXIT
   CASE ( 'STEQRV', 'steqrv' ); CALLST( 1 ) = .TRUE.
   CASE ( 'STEVXA', 'stebza' ); CALLST( 2 ) = .TRUE.
   CASE ( 'STEVXI', 'stebzi' ); CALLST( 3 ) = .TRUE.
   CASE ( 'STEVXV', 'stebzv' ); CALLST( 4 ) = .TRUE.
   CASE ( 'STEDCI', 'stedci' ); CALLST( 5 ) = .TRUE.
   CASE ( 'STEGRA', 'stegra' ); CALLST( 6 ) = .TRUE.
   CASE ( 'STEGRI', 'stegri' ); CALLST( 7 ) = .TRUE.
   CASE ( 'STEGRV', 'stegrv' ); CALLST( 8 ) = .TRUE.
   CASE DEFAULT; CALL HANDLER( 1, MACRO )
   END SELECT
END DO
!
END SUBROUTINE SSTEMCCALLST
