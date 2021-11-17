SUBROUTINE SSTEMCMATRIX( MAXN, NCASE, MACRO, T )
!
USE GSTEDEFINITIONS
USE GSTEINTERFACES1, ONLY : GETINTGR, PARSER
USE GSTEINTERFACES2, ONLY : PARSERLIST
USE SSTEDEFINITIONS
! 
!.. Scalar Arguments ..
CHARACTER( LEN=MACRO_NAME_LENGTH ) :: MACRO
INTEGER :: MAXN, NCASE
!
!.. Derived Data Type Argument ..
TYPE( T_LIST ), POINTER :: T
!
!==============================================================================!
!                                                                              !
!  Purpose:                                                                    !
!  =======                                                                     !
!                                                                              !
!  SSTEMCMATRIX deals with the macro that sets matrices to be generated.       !
!                                                                              !
!  Arguments:                                                                  !
!  =========                                                                   !
!                                                                              !
!  MAXN   (input/output) INTEGER                                               !
!         Maximum dimension of the matrix to be tested.                        !
!                                                                              !
!  NCASE  (input/output) INTEGER                                               !
!         Number of tridiagonal matrices defined in T.                         !
!                                                                              !
!  MACRO  (input) CHARACTER                                                    !
!         Macro definition.                                                    !
!                                                                              !
!  T      (input/output) T_LIST (derived data type)                            !
!         Properties of the tridiagonal matrices to be used in the tests.      !
!                                                                              !
!==============================================================================!
! 
!.. Local Scalars ..
CHARACTER( LEN=RECORD_LENGTH ) :: RECORD
INTEGER :: I, IERR, J, NMSIZE, NMTYPE
!
!.. Allocatable Arrays ..
INTEGER, ALLOCATABLE :: MSIZE( : ), MTYPE( : )
!
!.. Derived data types ..
TYPE( DATA_FROM_RECORD ), POINTER :: LIST
!
!.. External Subroutine ..
EXTERNAL HANDLER
!
!.. External Functions ..
CHARACTER( LEN=RECORD_LENGTH ), EXTERNAL :: GETRECORD
INTEGER, EXTERNAL :: LISTLENGHT
!
!.. Intrinsic Function ..
INTRINSIC MAX
!
!.. Executable Statements ......................................................
!
! Read matrix types.
!
RECORD = GETRECORD( )
LIST => PARSER( RECORD ); NMTYPE = LISTLENGHT( LIST )
ALLOCATE( MTYPE( NMTYPE ), STAT=IERR )
IF ( IERR /= 0 ) CALL HANDLER( 2, 'MTYPE (subroutine SSTEMCMATRIX)' )
MTYPE = GETINTGR( 'MTYPE', LIST, NMTYPE )
!
! Read matrix sizes.
!
RECORD = GETRECORD( )
LIST => PARSER( RECORD )
LIST => PARSERLIST( MACRO // '(invalid matrix size)', LIST )
NMSIZE = LISTLENGHT( LIST )
ALLOCATE( MSIZE( NMSIZE ), STAT=IERR )
IF ( IERR /= 0 ) CALL HANDLER( 2, 'MSIZE (subroutine SSTEMCMATRIX)' )
MSIZE = GETINTGR( 'MSIZE', LIST, NMSIZE )
!
! Check data consistency.
!
IF      ( MINVAL( MTYPE ) < 0 ) THEN
        CALL HANDLER( 1, MACRO // '(matrix type < 0)' )
ELSE IF ( MINVAL( MSIZE ) < 0 ) THEN
        CALL HANDLER( 1, MACRO // '(matrix size < 0' )
END IF
!
! Store information in data type T_DATA.
!
DO I = 1, NMTYPE
   DO J = 1, NMSIZE
      T%DATA = T_DATA( 2, MTYPE( I ), MSIZE( J ), 1, 1, 0, ISEED_INIT, ZERO )
      ALLOCATE( T%NEXT);  T => T%NEXT; NULLIFY( T%NEXT )
      MAXN = MAX( MAXN, MSIZE( J ) )
   END DO
END DO
!
NCASE = NCASE + NMTYPE*NMSIZE
!
! Deallocate arrays.
!
DEALLOCATE( MTYPE, MSIZE )
!
END SUBROUTINE SSTEMCMATRIX
