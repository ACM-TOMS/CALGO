      module mod_rmbv_csr
! **********************************************************************
!     Author : C. Voemel
!
!     Date of last modification : 7.7.00
!      
!     Description : PERFORMS MV MULT. WITH MATRIX IN 'CSR'-STORAGE
!                   rmbv = Right Multiplication By Vector: y=Ax
! **********************************************************************

      use representation_of_data
      use properties
      implicit none
      
      interface rmbv_csr
        module procedure irmbv_csr
        module procedure srmbv_csr
        module procedure drmbv_csr
        module procedure crmbv_csr
        module procedure zrmbv_csr
      end interface
      
      contains
! **********************************************************************
!
! **********************************************************************
#include "integerrmbv_csr_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singlePrecisionrmbv_csr_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doublePrecisionrmbv_csr_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singleComplexrmbv_csr_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doubleComplexrmbv_csr_source.F"
! **********************************************************************
!
! **********************************************************************
      end module mod_rmbv_csr
