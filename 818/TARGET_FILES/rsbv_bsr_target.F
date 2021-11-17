      module mod_rsbv_bsr
! **********************************************************************
!     Author : C. Voemel
!
!     Date of last modification : 7.7.00
!      
!     Description : PERFORMS TRI. SOLVE WITH MATRIX IN 'BSR'-STORAGE
!                   rsbv = Right Solve By Vector
! **********************************************************************
      use representation_of_data
      use properties
      use mod_dense_mat_algos

      implicit none

      interface rsbv_bsr
        module procedure irsbv_bsr
        module procedure srsbv_bsr
        module procedure drsbv_bsr
        module procedure crsbv_bsr
        module procedure zrsbv_bsr
      end interface

      contains
! **********************************************************************
!
! **********************************************************************
#include "integerrsbv_bsr_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singlePrecisionrsbv_bsr_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doublePrecisionrsbv_bsr_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singleComplexrsbv_bsr_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doubleComplexrsbv_bsr_source.F"
! **********************************************************************
!
! **********************************************************************
      end module mod_rsbv_bsr
