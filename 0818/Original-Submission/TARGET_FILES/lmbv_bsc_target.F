      module mod_lmbv_bsc
! **********************************************************************
!     Author : C. Voemel
!
!     Date of last modification : 7.7.00
!      
!     Description : PERFORMS MV MULT. WITH MATRIX IN 'BSC'-STORAGE
!                   lmbv = Left Multiplication By Vector: y^T=x^TA
! **********************************************************************

      use representation_of_data
      use properties
      use mod_dense_mat_algos

      implicit none
      
      interface lmbv_bsc
        module procedure ilmbv_bsc
        module procedure slmbv_bsc
        module procedure dlmbv_bsc
        module procedure clmbv_bsc
        module procedure zlmbv_bsc
      end interface
      
      contains
! **********************************************************************
!
! **********************************************************************
#include "integerlmbv_bsc_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singlePrecisionlmbv_bsc_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doublePrecisionlmbv_bsc_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singleComplexlmbv_bsc_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doubleComplexlmbv_bsc_source.F"
! **********************************************************************
!
! **********************************************************************
      end module mod_lmbv_bsc
