      module mod_lsbv_csr
! **********************************************************************
!     Author : C. Voemel
!
!     Date of last modification : 7.7.00
!      
!     Description : PERFORMS TRI. SOLVE WITH TRANSPOSE IN 'CSR'-STORAGE
!                   lsbv = Left Solve By Vector
! **********************************************************************      
      use representation_of_data
      use properties
      
      implicit none

      interface lsbv_csr
        module procedure ilsbv_csr
        module procedure slsbv_csr
        module procedure dlsbv_csr
        module procedure clsbv_csr
        module procedure zlsbv_csr
      end interface

      contains
! **********************************************************************
!
! **********************************************************************
#include "integerlsbv_csr_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singlePrecisionlsbv_csr_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doublePrecisionlsbv_csr_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singleComplexlsbv_csr_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doubleComplexlsbv_csr_source.F"
! **********************************************************************
!
! **********************************************************************
      end module mod_lsbv_csr
