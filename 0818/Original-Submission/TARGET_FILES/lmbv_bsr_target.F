      module mod_lmbv_bsr
! **********************************************************************
!     Author : C. Voemel
!
!     Date of last modification : 7.7.00
!      
!     Description : PERFORMS MV MULT. WITH MATRIX IN 'BSR'-STORAGE
!                   lmbv = Left Multiplication By Vector: y^T=x^TA
! **********************************************************************

      use representation_of_data
      use properties
      use mod_dense_mat_algos

      implicit none
      
      interface lmbv_bsr
        module procedure ilmbv_bsr
        module procedure slmbv_bsr
        module procedure dlmbv_bsr
        module procedure clmbv_bsr
        module procedure zlmbv_bsr
      end interface
      
      contains
! **********************************************************************
!
! **********************************************************************
#include "integerlmbv_bsr_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singlePrecisionlmbv_bsr_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doublePrecisionlmbv_bsr_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singleComplexlmbv_bsr_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doubleComplexlmbv_bsr_source.F"
! **********************************************************************
!
! **********************************************************************
      end module mod_lmbv_bsr
