MODULE SSTEINTERFACES3
!
!==============================================================================!
!                                                                              !
! This module defines interfaces for top level subroutines.                    !
!                                                                              !
!==============================================================================!
!
!-------------------------------------------------------------------------------
INTERFACE SSTEGETDATA
   SUBROUTINE SSTEGETDATA( HBANDA, HBANDR, MAXN, NCASE, NRILIU, NRVLVU, MINTT, &
                           TOLEBZ, TOLEGR, CALLST, DUMP, ILIU, NILIU, &
                           VLVU, NVLVU, EIGV, MTRX, T )
      USE GSTEDEFINITIONS
      USE GSTEINTERFACES1
      USE SSTEDEFINITIONS
      USE SSTEINTERFACES1
      USE SSTEINTERFACES2
      LOGICAL :: CALLST( 8 ), DUMP( 7 )
      INTEGER :: HBANDA, HBANDR, MAXN, NCASE, NILIU, NRILIU, NRVLVU, NVLVU
      REAL( KIND=PREC ) :: MINTT, TOLEBZ, TOLEGR
      TYPE( EIGV_LIST ), POINTER :: EIGV
      TYPE( MTRX_LIST ), POINTER :: MTRX
      TYPE( I_LIST ), POINTER :: ILIU
      TYPE( V_LIST ), POINTER :: VLVU
      TYPE( T_LIST ), POINTER :: T
   END SUBROUTINE SSTEGETDATA
END INTERFACE
!-------------------------------------------------------------------------------
INTERFACE SSTERUNTESTS
   SUBROUTINE SSTERUNTESTS( HBANDA, HBANDR, MAXN, NCASE, NRILIU, NRVLVU, MINTT, &
                            TOLEBZ, TOLEGR, CALLST, DUMP, ILIU, NILIU, &
                            VLVU, NVLVU, EIGV, MTRX, T )
      USE GSTEDEFINITIONS
      USE SSTEDEFINITIONS
      LOGICAL :: CALLST( 8 ), DUMP( 7 )
      INTEGER :: HBANDA, HBANDR, MAXN, NCASE, NILIU, NRILIU, NRVLVU, NVLVU
      REAL( KIND=PREC ) :: MINTT, TOLEBZ, TOLEGR
      TYPE( EIGV_LIST ), POINTER :: EIGV
      TYPE( MTRX_LIST ), POINTER :: MTRX
      TYPE( I_LIST ), POINTER :: ILIU
      TYPE( V_LIST ), POINTER :: VLVU
      TYPE( T_LIST ), POINTER :: T
   END SUBROUTINE SSTERUNTESTS
END INTERFACE
!-------------------------------------------------------------------------------
!
END MODULE SSTEINTERFACES3
