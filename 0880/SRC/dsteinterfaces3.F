MODULE DSTEINTERFACES3
!
!==============================================================================!
!                                                                              !
! This module defines interfaces for top level subroutines.                    !
!                                                                              !
!==============================================================================!
!
!-------------------------------------------------------------------------------
INTERFACE DSTEGETDATA
   SUBROUTINE DSTEGETDATA( HBANDA, HBANDR, MAXN, NCASE, NRILIU, NRVLVU, MINTT, &
                           TOLEBZ, TOLEGR, CALLST, DUMP, ILIU, NILIU, &
                           VLVU, NVLVU, EIGV, MTRX, T )
      USE GSTEDEFINITIONS
      USE GSTEINTERFACES1
      USE DSTEDEFINITIONS
      USE DSTEINTERFACES1
      USE DSTEINTERFACES2
      LOGICAL :: CALLST( 8 ), DUMP( 7 )
      INTEGER :: HBANDA, HBANDR, MAXN, NCASE, NILIU, NRILIU, NRVLVU, NVLVU
      REAL( KIND=PREC ) :: MINTT, TOLEBZ, TOLEGR
      TYPE( EIGV_LIST ), POINTER :: EIGV
      TYPE( MTRX_LIST ), POINTER :: MTRX
      TYPE( I_LIST ), POINTER :: ILIU
      TYPE( V_LIST ), POINTER :: VLVU
      TYPE( T_LIST ), POINTER :: T
   END SUBROUTINE DSTEGETDATA
END INTERFACE
!-------------------------------------------------------------------------------
INTERFACE DSTERUNTESTS
   SUBROUTINE DSTERUNTESTS( HBANDA, HBANDR, MAXN, NCASE, NRILIU, NRVLVU, MINTT, &
                            TOLEBZ, TOLEGR, CALLST, DUMP, ILIU, NILIU, &
                            VLVU, NVLVU, EIGV, MTRX, T )
      USE GSTEDEFINITIONS
      USE DSTEDEFINITIONS
      LOGICAL :: CALLST( 8 ), DUMP( 7 )
      INTEGER :: HBANDA, HBANDR, MAXN, NCASE, NILIU, NRILIU, NRVLVU, NVLVU
      REAL( KIND=PREC ) :: MINTT, TOLEBZ, TOLEGR
      TYPE( EIGV_LIST ), POINTER :: EIGV
      TYPE( MTRX_LIST ), POINTER :: MTRX
      TYPE( I_LIST ), POINTER :: ILIU
      TYPE( V_LIST ), POINTER :: VLVU
      TYPE( T_LIST ), POINTER :: T
   END SUBROUTINE DSTERUNTESTS
END INTERFACE
!-------------------------------------------------------------------------------
!
END MODULE DSTEINTERFACES3
