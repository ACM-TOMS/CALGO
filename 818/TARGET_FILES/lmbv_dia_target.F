      module mod_lmbv_dia
! **********************************************************************
!     Author : C. Voemel
!
!     Date of last modification : 7.7.00
!      
!     Description : PERFORMS MV MULT. WITH MATRIX IN 'DIA'-STORAGE
!                   lmbv = Left Multiplication By Vector: y^T=x^TA
! **********************************************************************

      use representation_of_data
      use properties
      implicit none
      
      interface lmbv_dia
        module procedure ilmbv_dia
        module procedure slmbv_dia
        module procedure dlmbv_dia
        module procedure clmbv_dia
        module procedure zlmbv_dia
      end interface
      
      contains
! **********************************************************************
!
! **********************************************************************
#include "integerlmbv_dia_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singlePrecisionlmbv_dia_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doublePrecisionlmbv_dia_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singleComplexlmbv_dia_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doubleComplexlmbv_dia_source.F"
! **********************************************************************
!
! **********************************************************************
      end module mod_lmbv_dia
