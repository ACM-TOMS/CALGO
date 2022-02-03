      module mod_lsbv_bco
! **********************************************************************
!     Author : C. Voemel
!
!     Date of last modification : 7.7.00
!      
!     Description : PERFORMS TRI. SOLVE WITH MATRIX IN 'BCO'-STORAGE
!                   lsbv = Left Solve By Vector
! **********************************************************************
      use mod_hash
      use representation_of_data
      use properties
      use mod_dense_mat_algos

      implicit none

      interface lsbv_bco
        module procedure ilsbv_bco
        module procedure slsbv_bco
        module procedure dlsbv_bco
        module procedure clsbv_bco
        module procedure zlsbv_bco
      end interface

      contains
! **********************************************************************
!
! **********************************************************************
#include "integerlsbv_bco_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singlePrecisionlsbv_bco_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doublePrecisionlsbv_bco_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singleComplexlsbv_bco_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doubleComplexlsbv_bco_source.F"
! **********************************************************************
!
! **********************************************************************
      end module mod_lsbv_bco
