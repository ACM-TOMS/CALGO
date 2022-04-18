      module mod_lmbv_bco
! **********************************************************************
!     Author : C. Voemel
!
!     Date of last modification : 7.7.00
!      
!     Description : PERFORMS MV MULT. WITH MATRIX IN 'BCO'-STORAGE
!                   lmbv = Left Multiplication By Vector: y^T=x^TA
! **********************************************************************

      use representation_of_data
      use properties
      use mod_dense_mat_algos

      implicit none
      
      interface lmbv_bco
        module procedure ilmbv_bco
        module procedure slmbv_bco
        module procedure dlmbv_bco
        module procedure clmbv_bco
        module procedure zlmbv_bco
      end interface
      
      contains
! **********************************************************************
!
! **********************************************************************
#include "integerlmbv_bco_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singlePrecisionlmbv_bco_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doublePrecisionlmbv_bco_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singleComplexlmbv_bco_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doubleComplexlmbv_bco_source.F"
! **********************************************************************
!
! **********************************************************************
      end module mod_lmbv_bco
