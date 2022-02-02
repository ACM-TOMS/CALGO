      module mod_rsbv_bsc
! **********************************************************************
!     Author : C. Voemel
!
!     Date of last modification : 7.7.00
!      
!     Description : PERFORMS TRI. SOLVE WITH MATRIX IN 'BSC'-STORAGE
!                   rsbv = Right Solve By Vector
! **********************************************************************
      use representation_of_data
      use properties
      use mod_dense_mat_algos

      implicit none

      interface rsbv_bsc
        module procedure irsbv_bsc
        module procedure srsbv_bsc
        module procedure drsbv_bsc
        module procedure crsbv_bsc
        module procedure zrsbv_bsc
      end interface

      contains
! **********************************************************************
!
! **********************************************************************
#include "integerrsbv_bsc_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singlePrecisionrsbv_bsc_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doublePrecisionrsbv_bsc_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singleComplexrsbv_bsc_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doubleComplexrsbv_bsc_source.F"
! **********************************************************************
!
! **********************************************************************
      end module mod_rsbv_bsc
