      module mod_lmbv_csr
! **********************************************************************
!     Author : C. Voemel
!
!     Date of last modification : 7.7.00
!      
!     Description : PERFORMS MV MULT. WITH MATRIX IN 'CSR'-STORAGE
!                   lmbv = Left Multiplication By Vector: y^T=x^TA
! **********************************************************************

      use representation_of_data
      use properties
      implicit none
      
      interface lmbv_csr
        module procedure ilmbv_csr
        module procedure slmbv_csr
        module procedure dlmbv_csr
        module procedure clmbv_csr
        module procedure zlmbv_csr
      end interface
      
      contains
! **********************************************************************
!
! **********************************************************************
#include "integerlmbv_csr_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singlePrecisionlmbv_csr_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doublePrecisionlmbv_csr_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singleComplexlmbv_csr_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doubleComplexlmbv_csr_source.F"
! **********************************************************************
!
! **********************************************************************
      end module mod_lmbv_csr
