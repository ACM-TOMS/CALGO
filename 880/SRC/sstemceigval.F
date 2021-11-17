SUBROUTINE SSTEMCEIGVAL( ECOND, EDIST, ESIGN, ISEED, MAXN, NCASE, MACRO, T )
!
USE GSTEDEFINITIONS
USE GSTEINTERFACES1, ONLY : GETINTGR, PARSER
USE GSTEINTERFACES2, ONLY : PARSERLIST
USE SSTEDEFINITIONS
! 
!.. Scalar Arguments ..
CHARACTER( LEN=MACRO_NAME_LENGTH ) :: MACRO
INTEGER :: ECOND, EDIST, ESIGN, MAXN, NCASE
!
!.. Array Argument ..
INTEGER :: ISEED( 4 )
!
!.. Derived Data Type Argument ..
TYPE( T_LIST ), POINTER :: T
!
!==============================================================================!
!                                                                              !
!  Purpose:                                                                    !
!  =======                                                                     !
!                                                                              !
!  SSTEMCEIGVAL deals with the macro that sets eigenvalue distributions.       !
!                                                                              !
!  Arguments:                                                                  !
!  =========                                                                   !
!                                                                              !
!  ECOND  (input) INTEGER                                                      !
!         Sets the condition number for ETYPE = 1:5,                           !
!         = 1, COND = 1 / SQRT(ULP), default                                   !
!         = 2, COND = 1 / (N*SQRT(ULP)),                                       !
!         = 3, COND = 1 / (10*N*SQRT(ULP)),                                    !
!         = 4, COND = 1 / ULP,                                                 !
!         = 5, COND = 1 / (N*ULP),                                             !
!         = 6, COND = 1 / (10*N*ULP),                                          !
!                                                                              !
!  EDIST  (input) INTEGER                                                      !
!         Specifies the type of the distribution to be used in the generation  !
!         of random eigenvalues                                                !
!         = 1, UNIFORM(  0, 1 )                                                !
!         = 2, UNIFORM( -1, 1 )                                                !
!         = 3, NORMAL( 0, 1 )                                                  !
!                                                                              !
!  ESIGN  (input) INTEGER                                                      !
!         Attributes signs to the entries in W                                 !
!         = 0, the entries of W are unchanged                                  !
!         = 1, multiplies each entry of W by +/- 1 with probability 0.5        !
!                                                                              !
!  ISEED  (input/output) INTEGER array, dimension ( 4 )                        !
!         Seed for the random number generator. Each entry of ISEED should     !
!         lie between 0 and 4095 inclusive and ISEED(4) should be odd.         !
!                                                                              !
!  MAXN   (input) INTEGER                                                      !
!         Maximum dimension of the matrix to be tested.                        !
!                                                                              !
!  NCASE  (output) INTEGER                                                     !
!         Number of tridiagonal matrices defined in T.                         !
!                                                                              !
!  MACRO  (input) CHARACTER                                                    !
!         Macro definition.                                                    !
!                                                                              !
!  T      (output) T_LIST (derived data type)                                  !
!         Properties of the tridiagonal matrices to be used in the tests.      !
!                                                                              !
!==============================================================================!
! 
!.. Local Scalars ..
CHARACTER( LEN=RECORD_LENGTH ) :: RECORD
INTEGER :: I, IERR, J, NESIZE, NETYPE
REAL( KIND=PREC ) :: TEMP
!
!.. Allocatable Arrays ..
INTEGER, ALLOCATABLE :: ESIZE( : ), ETYPE( : )
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
REAL( KIND=PREC ), EXTERNAL :: SLARAN
!
!.. Intrinsic Function ..
INTRINSIC MAX
!
!.. Executable Statements ......................................................
!
! Read distribution types.
!
RECORD = GETRECORD( )
LIST => PARSER( RECORD ); NETYPE = LISTLENGHT( LIST )
ALLOCATE( ETYPE( NETYPE ), STAT=IERR )
IF ( IERR /= 0 ) CALL HANDLER( 2, 'ETYPE (subroutine SSTEMCEIGVAL)' ) 
ETYPE = GETINTGR( MACRO, LIST, NETYPE )
!
! Read distribution sizes.
!
RECORD = GETRECORD( )
LIST => PARSER( RECORD )
LIST => PARSERLIST( MACRO // '(invalid matrix size)', LIST )
NESIZE = LISTLENGHT( LIST )
ALLOCATE( ESIZE( NESIZE ), STAT=IERR )
IF ( IERR /= 0 ) CALL HANDLER( 2, 'ESIZE (subroutine SSTEMCEIGVAL)' ) 
ESIZE = GETINTGR( MACRO, LIST, NESIZE )
!
! Check data consistency.
!
IF ( MINVAL( ESIZE ) < 0 ) THEN
   CALL HANDLER( 1, 'MATRIX (eigenvalue distribution size < 0)' )
END IF
!       
! Store information in data type T_LIST.
!
DO I = 1, NETYPE
   DO J = 1, NESIZE
      T%DATA = T_DATA( 1, ETYPE( I ), ESIZE( J ), ECOND, &
                       EDIST, ESIGN, ISEED, ZERO )
      ALLOCATE( T%NEXT);  T => T%NEXT; NULLIFY( T%NEXT )
      MAXN = MAX( MAXN, ESIZE( J ) )          
      TEMP = SLARAN( ISEED )
   END DO
END DO
!
NCASE = NCASE + NETYPE*NESIZE
!
! Deallocate arrays.
!
DEALLOCATE( ETYPE, ESIZE )
!
END SUBROUTINE SSTEMCEIGVAL
