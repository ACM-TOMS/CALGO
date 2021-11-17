      module mod_rsbv_bco
! **********************************************************************
!     Author : C. Voemel
!
!     Date of last modification : 7.7.00
!      
!     Description : PERFORMS TRI. SOLVE WITH MATRIX IN 'BCO'-STORAGE
!                   rsbv = Right Solve By Vector
! **********************************************************************
      use mod_hash
      use representation_of_data
      use properties
      use mod_dense_mat_algos

      implicit none

      interface rsbv_bco
        module procedure irsbv_bco
        module procedure srsbv_bco
        module procedure drsbv_bco
        module procedure crsbv_bco
        module procedure zrsbv_bco
      end interface

      contains
! **********************************************************************
!
! **********************************************************************
#include "integerrsbv_bco_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singlePrecisionrsbv_bco_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doublePrecisionrsbv_bco_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singleComplexrsbv_bco_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doubleComplexrsbv_bco_source.F"
! **********************************************************************
!
! **********************************************************************
      end module mod_rsbv_bco
