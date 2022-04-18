      module mod_lmbv_csc
! **********************************************************************
!     Author : C. Voemel
!
!     Date of last modification : 7.7.00
!      
!     Description : PERFORMS MV MULT. WITH MATRIX IN 'CSC'-STORAGE
!                   lmbv = Left Multiplication By Vector: y^T=x^TA
! **********************************************************************

      use representation_of_data
      use properties
      implicit none
      
      interface lmbv_csc
        module procedure ilmbv_csc
        module procedure slmbv_csc
        module procedure dlmbv_csc
        module procedure clmbv_csc
        module procedure zlmbv_csc
      end interface
      
      contains
! **********************************************************************
!
! **********************************************************************
#include "integerlmbv_csc_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singlePrecisionlmbv_csc_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doublePrecisionlmbv_csc_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singleComplexlmbv_csc_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doubleComplexlmbv_csc_source.F"
! **********************************************************************
!
! **********************************************************************
      end module mod_lmbv_csc
