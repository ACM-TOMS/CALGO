      module mod_rsbv_csr
! **********************************************************************
!     Author : C. Voemel
!
!     Date of last modification : 7.7.00
!      
!     Description : PERFORMS TRI. SOLVE WITH MATRIX IN 'CSR'-STORAGE
!                   rsbv = Right Solve By Vector
! **********************************************************************
      use representation_of_data
      use properties

      implicit none

      interface rsbv_csr
        module procedure irsbv_csr
        module procedure srsbv_csr
        module procedure drsbv_csr
        module procedure crsbv_csr
        module procedure zrsbv_csr
      end interface

      contains
! **********************************************************************
!
! **********************************************************************
#include "integerrsbv_csr_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singlePrecisionrsbv_csr_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doublePrecisionrsbv_csr_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singleComplexrsbv_csr_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doubleComplexrsbv_csr_source.F"
! **********************************************************************
!
! **********************************************************************
      end module mod_rsbv_csr
