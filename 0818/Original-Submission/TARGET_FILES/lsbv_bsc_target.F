      module mod_lsbv_bsc
! **********************************************************************
!     Author : C. Voemel
!
!     Date of last modification : 7.7.00
!      
!     Description : PERFORMS TRI. SOLVE WITH MATRIX IN 'BSC'-STORAGE
!                   lsbv = Left Solve By Vector
! **********************************************************************
      use representation_of_data
      use properties
      use mod_dense_mat_algos

      implicit none

      interface lsbv_bsc
        module procedure ilsbv_bsc
        module procedure slsbv_bsc
        module procedure dlsbv_bsc
        module procedure clsbv_bsc
        module procedure zlsbv_bsc
      end interface

      contains
! **********************************************************************
!
! **********************************************************************
#include "integerlsbv_bsc_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singlePrecisionlsbv_bsc_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doublePrecisionlsbv_bsc_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singleComplexlsbv_bsc_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doubleComplexlsbv_bsc_source.F"
! **********************************************************************
!
! **********************************************************************
      end module mod_lsbv_bsc
