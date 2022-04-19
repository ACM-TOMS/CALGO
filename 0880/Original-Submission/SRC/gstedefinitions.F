MODULE GSTEDEFINITIONS
!
! This module sets variables and defines type DATA_FROM_RECORD.
!
!*******************************************
!**                                       **
!**  The following lines may be modified  **
!**                                       **
!*******************************************
!
INTEGER, PARAMETER :: FUIN  = 5
INTEGER, PARAMETER :: FUOUT = 6
INTEGER, PARAMETER :: FUEXT = 15
INTEGER, PARAMETER :: FIELD_LENGTH      =  32
INTEGER, PARAMETER :: FILE_NAME_LENGTH  =  32 ! * >= FIELD_LENGTH *
INTEGER, PARAMETER :: MACRO_NAME_LENGTH =   8
INTEGER, PARAMETER :: RECORD_LENGTH     = 120
!
INTEGER, DIMENSION( 5 ) :: FUDUMP = (/ 16, 17, 18, 19, 20 /)
INTEGER, DIMENSION( 4 ) :: ISEED_INIT = (/ 4, 3, 2, 1 /)
!
CHARACTER( LEN=8 ) :: RFORMAT = '(E32.0) ' ! * set to FIELD_LENGTH *
CHARACTER( LEN=8 ) :: IFORMAT = '(I32)   ' ! * set to FIELD_LENGTH *
!
!*****************************************
!**                                     **
!**  Do not modify the following lines  **
!**                                     **
!*****************************************
!
CHARACTER( LEN=1 ), DIMENSION( 1:3 ) :: DIST = (/ 'U', 'S', 'N' /)
CHARACTER( LEN=1 ), DIMENSION( 0:1 ) :: SYMM = (/ 'P', 'S' /)
!
! List of data from record.
!
TYPE DATA_FROM_RECORD
     CHARACTER( LEN=FIELD_LENGTH ) :: FIELD
     TYPE( DATA_FROM_RECORD ), POINTER :: NEXT
END TYPE DATA_FROM_RECORD
!
END MODULE GSTEDEFINITIONS
