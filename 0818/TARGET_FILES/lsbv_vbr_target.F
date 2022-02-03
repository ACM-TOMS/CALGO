      module mod_lsbv_vbr
! **********************************************************************
!     Author : C. Voemel
!
!     Date of last modification : 7.7.00
!      
!     Description : PERFORMS TRI. SOLVE WITH MATRIX IN 'VBR'-STORAGE
!                   lsbv = Left Solve By Vector
! **********************************************************************
      use representation_of_data
      use properties
      use mod_dense_mat_algos

      implicit none

      interface lsbv_vbr
        module procedure ilsbv_vbr      
        module procedure slsbv_vbr
        module procedure dlsbv_vbr
        module procedure clsbv_vbr
        module procedure zlsbv_vbr
      end interface

      contains
! **********************************************************************
!
! **********************************************************************
#include "integerlsbv_vbr_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singlePrecisionlsbv_vbr_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doublePrecisionlsbv_vbr_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singleComplexlsbv_vbr_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doubleComplexlsbv_vbr_source.F"
! **********************************************************************
!
! **********************************************************************
      end module mod_lsbv_vbr
