SUBROUTINE SSTEMCEIGVV( MACRO, VLVU, NVLVU )
!
USE GSTEDEFINITIONS
USE GSTEINTERFACES1, ONLY : PARSER
USE SSTEDEFINITIONS
USE SSTEINTERFACES1, ONLY : GETSREAL
! 
!.. Scalar Arguments ..
CHARACTER( LEN=MACRO_NAME_LENGTH ) :: MACRO
INTEGER :: NVLVU
!
!.. Derived Data Type Argument ..
TYPE( V_LIST ), POINTER :: VLVU
!
!==============================================================================!
!                                                                              !
!  Purpose:                                                                    !
!  =======                                                                     !
!                                                                              !
!  SSTEMCEIGVV deals with the macro that sets intervals for the smallest and   !
!  largest eigenvalues to be computed.                                         !
!                                                                              !
!  Arguments:                                                                  !
!  =========                                                                   !
!                                                                              !
!  MACRO  (input) CHARACTER                                                    !
!         Macro definition.                                                    !
!                                                                              !
!  LIST   (input) DATA_FROM_RECORD (derived data type)                         !
!         String (file name).                                                  !
!                                                                              !
!  VLVU   (input/output) V_LIST (derived data type)                            !
!         Lower and upper bounds of intervals to be searched for eigenvalues.  !
!                                                                              !
!  NVLVU  (input/output) INTEGER                                               !
!         Number of intervals stored in VLVU.                                  !
!                                                                              !
!==============================================================================!
! 
!.. Local Scalars ..
CHARACTER( LEN=RECORD_LENGTH ) :: RECORD
INTEGER :: I, IERR, NVL, NVU
!
!.. Allocatable Arrays ..
REAL( KIND=PREC ), ALLOCATABLE :: VLLIST( : ), VULIST( : )
!
!.. Derived Data Type ..
TYPE( DATA_FROM_RECORD ), POINTER :: LIST
!
!.. External Subroutine ..
EXTERNAL HANDLER
!
!.. External Functions ..
CHARACTER( LEN=RECORD_LENGTH ), EXTERNAL :: GETRECORD
INTEGER, EXTERNAL :: LISTLENGHT
!
!.. Executable Statements ......................................................
!
RECORD = GETRECORD( )
LIST => PARSER( RECORD ); NVL = LISTLENGHT( LIST )
ALLOCATE( VLLIST( NVL ), STAT=IERR )
IF ( IERR /= 0 ) CALL HANDLER( 2, 'VLLIST (subroutine SSTEMCEIGVV)' ) 
VLLIST = GETSREAL( MACRO, LIST, NVL )
!
RECORD = GETRECORD( )
LIST => PARSER( RECORD ); NVU = LISTLENGHT( LIST )
ALLOCATE( VULIST( NVU ), STAT=IERR )
IF ( IERR /= 0 ) CALL HANDLER( 2, 'VULIST (subroutine SSTEMCEIGVV)' ) 
VULIST = GETSREAL( MACRO, LIST, NVU )
! 
IF ( NVL /= NVU ) CALL HANDLER( 1, MACRO // '(number of VL and VU differ)' )
!
DO I = 1, NVL
   IF ( VLLIST( I ) > VULIST( I ) ) THEN
      CALL HANDLER( 1, MACRO // '(VL > VU)' )
   ELSE
      VLVU%VL = VLLIST( I )
      VLVU%VU = VULIST( I )
      ALLOCATE( VLVU%NEXT);  VLVU => VLVU%NEXT; NULLIFY( VLVU%NEXT )
   END IF
END DO
!
NVLVU = NVLVU + NVL
!
END SUBROUTINE SSTEMCEIGVV
