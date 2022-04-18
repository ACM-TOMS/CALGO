      module mod_uscr_csr
! **********************************************************************
!     Author : C. Voemel
!
!     Date of last modification : 7.7.00
!      
!     Description : CREATION ROUTINE FOR MATRIX HANDLE FROM 'CSR'-FORMAT
! **********************************************************************

      use representation_of_data
      use properties
      use mod_usds

      implicit none
      
      interface uscr_csr
        module procedure iuscr_csr
        module procedure suscr_csr
        module procedure duscr_csr
        module procedure cuscr_csr
        module procedure zuscr_csr
      end interface
      
      contains
! **********************************************************************
!
! **********************************************************************
#include "integeruscr_csr_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singlePrecisionuscr_csr_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doublePrecisionuscr_csr_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singleComplexuscr_csr_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doubleComplexuscr_csr_source.F"
! **********************************************************************
!
! **********************************************************************
      end module mod_uscr_csr
