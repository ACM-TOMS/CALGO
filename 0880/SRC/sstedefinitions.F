MODULE SSTEDEFINITIONS
!
!==============================================================================!
!                                                                              !
! This module sets variables and defines derived data types.                   !
!                                                                              !
!==============================================================================!
!
!*******************************************
!**                                       **
!**  The following lines may be modified  **
!**                                       **
!*******************************************
!
INTEGER, PARAMETER :: PREC = KIND( 1.0E0 )
!
REAL( KIND=PREC ), PARAMETER :: HALF  = 0.5E0
REAL( KIND=PREC ), PARAMETER :: HNDRD = 1.0E2
REAL( KIND=PREC ), PARAMETER :: ONE   = 1.0E0
REAL( KIND=PREC ), PARAMETER :: TEN   = 1.0E1
REAL( KIND=PREC ), PARAMETER :: THSND = 1.0E3
REAL( KIND=PREC ), PARAMETER :: TWO   = 2.0E0
REAL( KIND=PREC ), PARAMETER :: ZERO  = 0.0E0
!
CHARACTER( LEN=18 ), DIMENSION( 1:5 ) :: FDUMP  = (/ 'sstetester.out.T  ', &
                                                     'sstetester.out.W  ', &
                                                     'sstetester.out.Z  ', &
                                                     'sstetester.out.log', &
                                                     'sstetester.out.m  ' /)
CHARACTER( LEN=15 ), DIMENSION( 1:8 ) :: IDTEST = (/ 'SSTEQR(COMPZ=I)', &
                                                     'SSTEVX(RANGE=A)', &
                                                     'SSTEVX(RANGE=I)', &
                                                     'SSTEVX(RANGE=V)', &
                                                     'SSTEDC(COMPZ=I)', &
                                                     'SSTEGR(RANGE=A)', &
                                                     'SSTEGR(RANGE=I)', &
                                                     'SSTEGR(RANGE=V)' /)
!
!*****************************************
!**                                     **
!**  Do not modify the following lines  **
!**                                     **
!*****************************************
!
! Derived data type for tridiagonals.
!
TYPE T_DATA
     INTEGER :: FORM, TYPE, SIZE, COND, DIST, SIGN, SEED( 4 )
     REAL( KIND=PREC ) :: EN
END TYPE T_DATA
!
! List of tridiagonals.
!
TYPE T_LIST
     TYPE( T_DATA ) :: DATA
     TYPE( T_LIST ), POINTER :: NEXT
END TYPE T_LIST
!
! List of eigenvalue distributions read from files.
!
TYPE EIGV_LIST
     REAL( KIND=PREC ), POINTER :: S( : )
     TYPE( EIGV_LIST ), POINTER :: NEXT
END TYPE EIGV_LIST
!
! List of tridiagonal matrices read from files.
!
TYPE MTRX_LIST
     REAL( KIND=PREC ), POINTER :: D( : ), E( : )
     TYPE( MTRX_LIST ), POINTER :: NEXT
END TYPE MTRX_LIST
!
! List of ranges, indexes.
!
TYPE I_LIST
     INTEGER :: IL, IU
     TYPE( I_LIST ), POINTER :: NEXT
END TYPE I_LIST
!
! List of ranges, values.
!
TYPE V_LIST
     REAL( KIND=PREC ) :: VL, VU
     TYPE( V_LIST ), POINTER :: NEXT
END TYPE V_LIST
!
END MODULE SSTEDEFINITIONS
