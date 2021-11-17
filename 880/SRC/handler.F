SUBROUTINE HANDLER( INFO, STRING )
!
USE GSTEDEFINITIONS
!
!.. Scalar Arguments ..
CHARACTER( * ) :: STRING
INTEGER :: INFO
!
!==============================================================================!
!                                                                              !
!  Purpose:                                                                    !
!  =======                                                                     !
!                                                                              !
!  HANDLER prints error messages and terminates execution.                     !
!                                                                              !
!  Arguments:                                                                  !
!  =========                                                                   !
!                                                                              !
!  INFO    (input) INTEGER                                                     !
!          Sets the error message to be printed, see cases below.              !
!                                                                              !
!  STRING  (input) CHARACTER                                                   !
!          Additional string to be printed.                                    !
!                                                                              !
!==============================================================================!
! 
!.. Local Scalar ..
CHARACTER( LEN=RECORD_LENGTH ) :: PHRASE
!
!.. Executable Statements ......................................................
!
SELECT CASE ( INFO )
CASE ( 1 )
     PHRASE = 'Invalid data in input file: ' // STRING
CASE ( 2 )
     PHRASE = 'Failure in array allocation: ' // STRING
CASE ( 3 )
     PHRASE = 'Parsing error in key word: ' // STRING
CASE ( 4 )
     PHRASE = 'Input file not found: ' // STRING
CASE ( 5 )
     PHRASE = 'Key word not implemented: ' // STRING
CASE ( 6 )
     PHRASE = 'Nothing to be tested: ' // STRING
CASE DEFAULT
     IF ( INFO < 0 ) THEN
        WRITE ( PHRASE, FMT='(2A,I5)' ) STRING, ' fatal error: INFO = ', INFO
     ELSE
        PHRASE = STRING
     END IF
END SELECT
!
WRITE ( UNIT=FUOUT, FMT='(A3,A60)' ) '-> ', PHRASE
STOP '* Execution terminated by handler *'
!
END SUBROUTINE HANDLER
