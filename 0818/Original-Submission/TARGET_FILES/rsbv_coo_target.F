      module mod_rsbv_coo
! **********************************************************************
!     Author : C. Voemel
!
!     Date of last modification : 7.7.00
!      
!     Description : PERFORMS TRI. SOLVE WITH MATRIX IN 'COO'-STORAGE
!                   rsbv = Right Solve By Vector
! **********************************************************************
      use mod_hash
      use representation_of_data
      use properties

      implicit none

      interface rsbv_coo
        module procedure irsbv_coo
        module procedure srsbv_coo
        module procedure drsbv_coo
        module procedure crsbv_coo
        module procedure zrsbv_coo
      end interface

      contains
! **********************************************************************
!
! **********************************************************************
#include "integerrsbv_coo_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singlePrecisionrsbv_coo_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doublePrecisionrsbv_coo_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singleComplexrsbv_coo_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doubleComplexrsbv_coo_source.F"
! **********************************************************************
!
! **********************************************************************
      end module mod_rsbv_coo
