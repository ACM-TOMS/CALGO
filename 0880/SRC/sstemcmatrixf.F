SUBROUTINE SSTEMCMATRIXF( MAXN, NCASE, MACRO, LIST, MTRX, T )
!
USE GSTEDEFINITIONS
USE GSTEINTERFACES1, ONLY : GETSTRNG
USE SSTEDEFINITIONS
! 
!.. Scalar Arguments ..
CHARACTER( LEN=MACRO_NAME_LENGTH ) :: MACRO
INTEGER :: MAXN, NCASE
!
!.. Derived Data Type Arguments ..
TYPE( DATA_FROM_RECORD ), POINTER :: LIST
TYPE( MTRX_LIST ), POINTER :: MTRX
TYPE( T_LIST ), POINTER :: T
!
!==============================================================================!
!                                                                              !
!  Purpose:                                                                    !
!  =======                                                                     !
!                                                                              !
!  SSTEMCMATRIXF deals with the macro that reads a matrix from a file.         !
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
!  LIST   (input) DATA_FROM_RECORD (derived data type)                         !
!         String (file name).                                                  !
!                                                                              !
!  MTRX   (input/output) MTRX_LIST (derived data type)                         !
!         Tridiagonal matrices read from files.                                !
!                                                                              !
!  T      (input/output) T_LIST (derived data type)                            !
!         Properties of the tridiagonal matrices to be used in the tests.      !
!                                                                              !
!==============================================================================!
! 
!.. Local Scalars ..
LOGICAL :: FILE_EXIST
INTEGER :: I, IERR, J, N
REAL( KIND=PREC ) :: TEMP1, TEMP2
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
   ALLOCATE( MTRX%D( N ), MTRX%E( N ), STAT=IERR )
   DO I = 1, N
      READ ( UNIT=FUEXT, FMT=*, IOSTAT=IERR ) J, TEMP1, TEMP2
      IF ( IERR /= 0 ) CALL HANDLER( 1, FILE_NAME( 1 ) )
      MTRX%D( J ) = TEMP1; MTRX%E( J ) = TEMP2
   END DO
   CLOSE( UNIT=FUEXT )
   T%DATA = T_DATA( 4, 1, N, 1, 1, 0, ISEED_INIT, ZERO )
   ALLOCATE( T%NEXT);  T => T%NEXT; NULLIFY( T%NEXT )
   ALLOCATE( MTRX%NEXT );  MTRX => MTRX%NEXT; NULLIFY( MTRX%NEXT )
ELSE
   CALL HANDLER( 4, FILE_NAME( 1 ) )
END IF
!
MAXN = MAX( MAXN, N )
NCASE = NCASE + 1
!
END SUBROUTINE SSTEMCMATRIXF
