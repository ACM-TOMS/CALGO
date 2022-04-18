      module mod_rsbv_vbr
! **********************************************************************
!     Author : C. Voemel
!
!     Date of last modification : 7.7.00
!      
!     Description : PERFORMS TRI. SOLVE WITH MATRIX IN 'VBR'-STORAGE
!                   rsbv = Right Solve By Vector
! **********************************************************************
      use representation_of_data
      use properties
      use mod_dense_mat_algos

      implicit none

      interface rsbv_vbr
        module procedure irsbv_vbr
        module procedure srsbv_vbr
        module procedure drsbv_vbr
        module procedure crsbv_vbr
        module procedure zrsbv_vbr
      end interface

      contains
! **********************************************************************
!
! **********************************************************************
#include "integerrsbv_vbr_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singlePrecisionrsbv_vbr_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doublePrecisionrsbv_vbr_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singleComplexrsbv_vbr_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doubleComplexrsbv_vbr_source.F"
! **********************************************************************
!
! **********************************************************************
      end module mod_rsbv_vbr
