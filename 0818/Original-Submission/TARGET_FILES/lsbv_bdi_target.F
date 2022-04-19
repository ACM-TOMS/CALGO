      module mod_lsbv_bdi
! **********************************************************************
!     Author : C. Voemel
!
!     Date of last modification : 7.7.00
!      
!     Description : PERFORMS TRI. SOLVE WITH MATRIX IN 'BDI'-STORAGE
!                   lsbv = Left Solve By Vector
! **********************************************************************
      use representation_of_data
      use properties
      use mod_dense_mat_algos

      implicit none

      interface lsbv_bdi
        module procedure ilsbv_bdi
        module procedure slsbv_bdi
        module procedure dlsbv_bdi
        module procedure clsbv_bdi
        module procedure zlsbv_bdi
      end interface

      contains
! **********************************************************************
!
! **********************************************************************
#include "integerlsbv_bdi_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singlePrecisionlsbv_bdi_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doublePrecisionlsbv_bdi_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singleComplexlsbv_bdi_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doubleComplexlsbv_bdi_source.F"
! **********************************************************************
!
! **********************************************************************
      end module mod_lsbv_bdi
