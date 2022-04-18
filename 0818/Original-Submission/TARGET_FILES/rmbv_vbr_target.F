      module mod_rmbv_vbr
! **********************************************************************
!     Author : C. Voemel
!
!     Date of last modification : 7.7.00
!      
!     Description : PERFORMS MV MULT. WITH MATRIX IN 'VBR'-STORAGE
!                   rmbv = Right Multiplication By Vector: y=Ax
! **********************************************************************

      use representation_of_data
      use properties
      use mod_dense_mat_algos

      implicit none
      
      interface rmbv_vbr
        module procedure irmbv_vbr
        module procedure srmbv_vbr
        module procedure drmbv_vbr
        module procedure crmbv_vbr
        module procedure zrmbv_vbr
      end interface
      
      contains
! **********************************************************************
!
! **********************************************************************
#include "integerrmbv_vbr_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singlePrecisionrmbv_vbr_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doublePrecisionrmbv_vbr_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singleComplexrmbv_vbr_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doubleComplexrmbv_vbr_source.F"
! **********************************************************************
!
! **********************************************************************
      end module mod_rmbv_vbr
