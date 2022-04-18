      module mod_lmbv_vbr
! **********************************************************************
!     Author : C. Voemel
!
!     Date of last modification : 7.7.00
!      
!     Description : PERFORMS MV MULT. WITH MATRIX IN 'VBR'-STORAGE
!                   lmbv = Left Multiplication By Vector: y^T=x^TA
! **********************************************************************

      use representation_of_data
      use properties
      use mod_dense_mat_algos

      implicit none
      
      interface lmbv_vbr
        module procedure ilmbv_vbr
        module procedure slmbv_vbr
        module procedure dlmbv_vbr
        module procedure clmbv_vbr
        module procedure zlmbv_vbr
      end interface
      
      contains
! **********************************************************************
!
! **********************************************************************
#include "integerlmbv_vbr_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singlePrecisionlmbv_vbr_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doublePrecisionlmbv_vbr_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singleComplexlmbv_vbr_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doubleComplexlmbv_vbr_source.F"
! **********************************************************************
!
! **********************************************************************
      end module mod_lmbv_vbr
