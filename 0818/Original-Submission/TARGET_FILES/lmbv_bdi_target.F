      module mod_lmbv_bdi
! **********************************************************************
!     Author : C. Voemel
!
!     Date of last modification : 7.7.00
!      
!     Description : PERFORMS MV MULT. WITH MATRIX IN 'BDI'-STORAGE
!                   lmbv = Left Multiplication By Vector: y^T=x^TA
! **********************************************************************

      use representation_of_data
      use properties
      use mod_dense_mat_algos

      implicit none
      
      interface lmbv_bdi
        module procedure ilmbv_bdi
        module procedure slmbv_bdi
        module procedure dlmbv_bdi
        module procedure clmbv_bdi
        module procedure zlmbv_bdi
      end interface
      
      contains
! **********************************************************************
!
! **********************************************************************
#include "integerlmbv_bdi_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singlePrecisionlmbv_bdi_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doublePrecisionlmbv_bdi_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singleComplexlmbv_bdi_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doubleComplexlmbv_bdi_source.F"
! **********************************************************************
!
! **********************************************************************
      end module mod_lmbv_bdi
