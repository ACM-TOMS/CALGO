      module mod_lmbv_coo
! **********************************************************************
!     Author : C. Voemel
!
!     Date of last modification : 7.7.00
!      
!     Description : PERFORMS MV MULT. WITH TRANSPOSE IN 'COO'-STORAGE
!                   lmbv = Left Multiplication By Vector: y^T = x^TA
! **********************************************************************

      use representation_of_data
      use properties
      implicit none
      
      interface lmbv_coo
        module procedure ilmbv_coo
        module procedure slmbv_coo
        module procedure dlmbv_coo
        module procedure clmbv_coo
        module procedure zlmbv_coo
      end interface
      
      contains
! **********************************************************************
!
! **********************************************************************
#include "integerlmbv_coo_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singlePrecisionlmbv_coo_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doublePrecisionlmbv_coo_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singleComplexlmbv_coo_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doubleComplexlmbv_coo_source.F"
! **********************************************************************
!
! **********************************************************************
      end module mod_lmbv_coo
