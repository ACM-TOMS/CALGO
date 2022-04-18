SUBROUTINE DSTEMCEIGVI( MACRO, ILIU, NILIU )
!
USE GSTEDEFINITIONS
USE GSTEINTERFACES1, ONLY : GETINTGR, PARSER
USE DSTEDEFINITIONS
! 
!.. Scalar Arguments ..
CHARACTER( LEN=MACRO_NAME_LENGTH ) :: MACRO
INTEGER :: NILIU
!
!.. Derived Data Type Argument ..
TYPE( I_LIST ), POINTER :: ILIU
!
!==============================================================================!
!                                                                              !
!  Purpose:                                                                    !
!  =======                                                                     !
!                                                                              !
!  DSTEMCEIGVI deals with the macro that sets indices of the smallest and      !
!  largest eigenvalues to be computed.                                         !
!                                                                              !
!  Arguments:                                                                  !
!  =========                                                                   !
!                                                                              !
!  MACRO  (input) CHARACTER                                                    !
!         Macro definition.                                                    !
!                                                                              !
!  ILIU   (input/output) I_LIST (derived data type)                            !
!         Indices (in ascending order) of the smallest and largest             !
!         eigenvalues to be computed.                                          !
!                                                                              !
!  NILIU  (input/output) INTEGER                                               !
!         Number of pair of indices stored in ILIU.                            !
!                                                                              !
!==============================================================================!
! 
!.. Local Scalars ..
CHARACTER( LEN=RECORD_LENGTH ) :: RECORD
INTEGER :: I, IERR, NIL, NIU
!
!.. Allocatable Arrays ..
INTEGER, ALLOCATABLE :: ILLIST( : ), IULIST( : )
!
!.. Derived Data Types ..
TYPE( DATA_FROM_RECORD ), POINTER :: LIST
!
!.. External Functions ..
CHARACTER( LEN=RECORD_LENGTH ), EXTERNAL :: GETRECORD
INTEGER, EXTERNAL :: LISTLENGHT
!
!.. Intrinsic Function ..
INTRINSIC MIN
!
!.. Executable Statements ......................................................
!
RECORD = GETRECORD( )
LIST => PARSER( RECORD ); NIL = LISTLENGHT( LIST )
ALLOCATE( ILLIST( NIL ), STAT=IERR )
IF ( IERR /= 0 ) CALL HANDLER( 2, 'ILLIST (subroutine DSTEMCEIGVI)' ) 
ILLIST = GETINTGR( MACRO, LIST, NIL )
!
RECORD = GETRECORD( )
LIST => PARSER( RECORD ); NIU = LISTLENGHT( LIST )
ALLOCATE( IULIST( NIU ), STAT=IERR )
IF ( IERR /= 0 ) CALL HANDLER( 2, 'IULIST (subroutine DSTEMCEIGVI)' ) 
IULIST = GETINTGR( MACRO, LIST, NIU )
! 
IF ( NIL /= NIU ) CALL HANDLER( 1, MACRO // '(number of IL and IU differ)' )
!
DO I = 1, NIL
   IF      ( MIN( ILLIST( I ), IULIST( I ) ) < 1 ) THEN
           CALL HANDLER( 1, MACRO // '(min(IL,IU) < 0)' )
   ELSE IF ( ILLIST( I ) > IULIST( I ) ) THEN
           CALL HANDLER( 1, MACRO // '(IL > IU)' )
   ELSE
           ILIU%IL = ILLIST( I )
           ILIU%IU = IULIST( I )
           ALLOCATE( ILIU%NEXT );  ILIU => ILIU%NEXT; NULLIFY( ILIU%NEXT )
   END IF
END DO
!
NILIU = NILIU + NIL
!
END SUBROUTINE DSTEMCEIGVI
