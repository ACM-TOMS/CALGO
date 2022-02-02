SUBROUTINE DSTEMCDUMP( MACRO, LIST, DUMP )
!
USE GSTEDEFINITIONS
USE GSTEINTERFACES1, ONLY : GETSTRNG
! 
!.. Scalar Argument ..
CHARACTER( LEN=MACRO_NAME_LENGTH ) :: MACRO
!
!.. Array Argument ..
LOGICAL :: DUMP( 7 )
!
!.. Derived Data Type Argument ..
TYPE( DATA_FROM_RECORD ), POINTER :: LIST
!
!==============================================================================!
!                                                                              !
!  Purpose:                                                                    !
!  =======                                                                     !
!                                                                              !
!  DSTEMCDUMP deals with the macro that sets variables to be saved.            !
!                                                                              !
!  Arguments:                                                                  !
!  =========                                                                   !
!                                                                              !
!  MACRO  (input) CHARACTER                                                    !
!         Macro definition.                                                    !
!                                                                              !
!  LIST   (input) DATA_FROM_RECORD (derived data type)                         !
!         List of strings.                                                     !
!                                                                              !
!  DUMP   (output) LOGICAL array, dimension ( 7 )                              !
!         Defines data to be written into files,                               !
!         DUMP( 1 ) : tridiagonal matrix (i,d_i,e_i)                           !
!         DUMP( 2 ) : eigenvalues                                              !
!         DUMP( 3 ) : eigenvectors                                             !
!         DUMP( 4 ) : timing, residuals, orthogonality                         !
!         DUMP( 5 ) : tridiagonal matrix (i,d_i,e_i) in Matlab format          !
!         DUMP( 6 ) : eigenvalues in Matlab format                             !
!         DUMP( 7 ) : eigenvectors in Matlab format                            !
!                                                                              !
!==============================================================================!
! 
!.. Local Scalars ..
INTEGER :: I, NDUMP
!
!.. Allocatable Array ..
CHARACTER( LEN=MACRO_NAME_LENGTH ), ALLOCATABLE :: LDUMP( : )
!
!.. External Subroutine ..
EXTERNAL HANDLER
!
!.. External Function ..
INTEGER, EXTERNAL :: LISTLENGHT
!
!.. Executable Statements ......................................................
!
NDUMP = LISTLENGHT( LIST )
ALLOCATE( LDUMP( NDUMP ) )
LDUMP = GETSTRNG( MACRO, LIST, NDUMP )
DO I = 1, NDUMP
   SELECT CASE ( LDUMP( I ) )
   CASE ( 'T'   , 't'   ); DUMP( 1 ) = .TRUE. 
   CASE ( 'W'   , 'w'   ); DUMP( 2 ) = .TRUE.
   CASE ( 'Z'   , 'z'   ); DUMP( 3 ) = .TRUE.
   CASE ( 'LOG' , 'log' ); DUMP( 4 ) = .TRUE.
   CASE ( 'T.M' , 't.m' ); DUMP( 5 ) = .TRUE.
   CASE ( 'W.M' , 'w.m' ); DUMP( 6 ) = .TRUE.
   CASE ( 'Z.M' , 'z.m' ); DUMP( 7 ) = .TRUE.
   CASE DEFAULT; CALL HANDLER( 1, MACRO )
   END SELECT
END DO
DEALLOCATE( LDUMP )
!
END SUBROUTINE DSTEMCDUMP
