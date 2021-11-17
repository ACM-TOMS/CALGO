      module mod_rsbv_bdi
! **********************************************************************
!     Author : C. Voemel
!
!     Date of last modification : 7.7.00
!      
!     Description : PERFORMS TRI. SOLVE WITH MATRIX IN 'BDI'-STORAGE
!                   rsbv = Right Solve By Vector
! **********************************************************************
      use representation_of_data
      use properties
      use mod_dense_mat_algos

      implicit none

      interface rsbv_bdi
        module procedure irsbv_bdi
        module procedure srsbv_bdi
        module procedure drsbv_bdi
        module procedure crsbv_bdi
        module procedure zrsbv_bdi
      end interface

      contains
! **********************************************************************
!
! **********************************************************************
#include "integerrsbv_bdi_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singlePrecisionrsbv_bdi_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doublePrecisionrsbv_bdi_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singleComplexrsbv_bdi_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doubleComplexrsbv_bdi_source.F"
! **********************************************************************
!
! **********************************************************************
      end module mod_rsbv_bdi
