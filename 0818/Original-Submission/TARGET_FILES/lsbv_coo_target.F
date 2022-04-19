      module mod_lsbv_coo
! **********************************************************************
!     Author : C. Voemel
!
!     Date of last modification : 7.7.00
!      
!     Description : PERFORMS TRI. SOLVE WITH MATRIX IN 'COO'-STORAGE
!                   lsbv = Left Solve By Vector
! **********************************************************************
      use mod_hash
      use representation_of_data
      use properties

      implicit none

      interface lsbv_coo
        module procedure ilsbv_coo
        module procedure slsbv_coo
        module procedure dlsbv_coo
        module procedure clsbv_coo
        module procedure zlsbv_coo
      end interface

      contains
! **********************************************************************
!
! **********************************************************************
#include "integerlsbv_coo_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singlePrecisionlsbv_coo_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doublePrecisionlsbv_coo_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singleComplexlsbv_coo_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doubleComplexlsbv_coo_source.F"
! **********************************************************************
!
! **********************************************************************
      end module mod_lsbv_coo
