SUBROUTINE DSTEMCEIGVALF( ISEED, MAXN, NCASE, MACRO, LIST, EIGV, T )
!
USE GSTEDEFINITIONS
USE GSTEINTERFACES1, ONLY : GETSTRNG
USE DSTEDEFINITIONS
! 
!.. Scalar Arguments ..
CHARACTER( LEN=MACRO_NAME_LENGTH ) :: MACRO
INTEGER :: MAXN, NCASE
!
!.. Array Argument ..
INTEGER :: ISEED( 4 )
!
!.. Derived Data Type Arguments ..
TYPE( DATA_FROM_RECORD ), POINTER :: LIST
TYPE( EIGV_LIST ), POINTER :: EIGV
TYPE( T_LIST ), POINTER :: T
!
!==============================================================================!
!                                                                              !
!  Purpose:                                                                    !
!  =======                                                                     !
!                                                                              !
!  DSTEMCEIGVALF deals with the macro that reads an eigenvalue distribution    !
!  from a file.                                                                !
!                                                                              !
!  Arguments:                                                                  !
!  =========                                                                   !
!                                                                              !
!  ISEED  (input/output) INTEGER array, dimension ( 4 )                        !
!         Seed for the random number generator. Each entry of ISEED should     !
!         lie between 0 and 4095 inclusive and ISEED(4) should be odd.         !
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
!  LIST   (input) DATA_FROM_RECORD (derived data type)                         !
!         String (file name).                                                  !
!                                                                              !
!  EIGV   (input/output) EIGV_LIST (derived data type)                         !
!         Eigenvalue distributions read from files.                            !
!                                                                              !
!  T      (input/output) T_LIST (derived data type)                            !
!         Properties of the tridiagonal matrices to be used in the tests.      !
!                                                                              !
!==============================================================================!
! 
!.. Local Scalars ..
LOGICAL :: FILE_EXIST
INTEGER :: I, IERR, N
REAL( KIND=PREC ) :: TEMP
! 
!.. Static Array ..
CHARACTER( LEN=FILE_NAME_LENGTH ) :: FILE_NAME( 1 )
!
!.. External Subroutine ..
EXTERNAL HANDLER
!
!.. Intrinsic Function ..
INTRINSIC MAX
!
!.. Executable Statements ......................................................
!
FILE_NAME = GETSTRNG( MACRO, LIST, 1 )
INQUIRE ( FILE=FILE_NAME( 1 ), EXIST=FILE_EXIST )
!
IF ( FILE_EXIST ) THEN
   OPEN( UNIT=FUEXT, FILE=FILE_NAME( 1 ), IOSTAT=IERR )
   IF ( IERR /= 0 ) CALL HANDLER( 4, FILE_NAME( 1 ) )
   READ ( UNIT=FUEXT, FMT=*, IOSTAT=IERR ) N
   IF ( IERR /= 0 ) CALL HANDLER( 1, FILE_NAME( 1 ) )
   ALLOCATE( EIGV%S( N ) )
   DO I = 1, N
      READ ( UNIT=FUEXT, FMT=*, IOSTAT=IERR ) TEMP
      IF ( IERR /= 0 ) CALL HANDLER( 1, FILE_NAME( 1 ) )
      EIGV%S( I ) = TEMP
   END DO
   CLOSE( UNIT=FUEXT )
   T%DATA = T_DATA( 3, 1, N, 1, 1, 0, ISEED, ZERO )
   ALLOCATE( T%NEXT );  T => T%NEXT; NULLIFY( T%NEXT )
   ALLOCATE( EIGV%NEXT );  EIGV => EIGV%NEXT; NULLIFY( EIGV%NEXT )
ELSE
   CALL HANDLER( 4, FILE_NAME( 1 ) )
END IF
!
MAXN = MAX( MAXN, N )
NCASE = NCASE + 1
!
END SUBROUTINE DSTEMCEIGVALF
