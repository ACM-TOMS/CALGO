RECURSIVE FUNCTION GETRECORD( ) RESULT( RECORD )
!
USE GSTEDEFINITIONS
!
!.. Function Result ..
CHARACTER( LEN=RECORD_LENGTH ) :: RECORD
!
!==============================================================================!
!                                                                              !
!  Purpose:                                                                    !
!  =======                                                                     !
!                                                                              !
!  GETRECORD reads a record from the input file.                               !
!                                                                              !
!==============================================================================!
!  
!.. Local Scalar ..
INTEGER :: IOERR
! 
!.. External Subroutine ..
EXTERNAL HANDLER
!
!.. Intrinsic Function ..
INTRINSIC ADJUSTL 
!
!.. Executable Statements ......................................................
!
! Read one record from the input file.
!
READ ( UNIT=FUIN, IOSTAT=IOERR, FMT='(A)' ) RECORD
!
IF      ( IOERR < 0 ) THEN   !** End of file **
        RECORD = 'END'
ELSE IF ( IOERR > 0 ) THEN   !** IO error, quit **
        CALL HANDLER( 0, 'IO error in input file' )
END IF 
!
! Adjust left (remove leading blanks).
!
RECORD = ADJUSTL( RECORD )
!
! If this is a comment line, read another record.     
!
IF ( RECORD( 1:1 ) == ' ' ) CALL HANDLER( 0, 'Empty line in input file' )
IF ( RECORD( 1:1 ) == '%' ) RECORD = GETRECORD( )
!
END FUNCTION GETRECORD
