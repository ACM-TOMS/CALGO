MODULE DSTEDEFINITIONS
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
INTEGER, PARAMETER :: PREC = KIND( 1.0D0 )
!
REAL( KIND=PREC ), PARAMETER :: HALF  = 0.5D0
REAL( KIND=PREC ), PARAMETER :: HNDRD = 1.0D2
REAL( KIND=PREC ), PARAMETER :: ONE   = 1.0D0
REAL( KIND=PREC ), PARAMETER :: TEN   = 1.0D1
REAL( KIND=PREC ), PARAMETER :: THSND = 1.0D3
REAL( KIND=PREC ), PARAMETER :: TWO   = 2.0D0
REAL( KIND=PREC ), PARAMETER :: ZERO  = 0.0D0
!
CHARACTER( LEN=18 ), DIMENSION( 1:5 ) :: FDUMP  = (/ 'dstetester.out.T  ', &
                                                     'dstetester.out.W  ', &
                                                     'dstetester.out.Z  ', &
                                                     'dstetester.out.log', &
                                                     'dstetester.out.m  ' /)
CHARACTER( LEN=15 ), DIMENSION( 1:8 ) :: IDTEST = (/ 'DSTEQR(COMPZ=I)', &
                                                     'DSTEVX(RANGE=A)', &
                                                     'DSTEVX(RANGE=I)', &
                                                     'DSTEVX(RANGE=V)', &
                                                     'DSTEDC(COMPZ=I)', &
                                                     'DSTEGR(RANGE=A)', &
                                                     'DSTEGR(RANGE=I)', &
                                                     'DSTEGR(RANGE=V)' /)
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
END MODULE DSTEDEFINITIONS
