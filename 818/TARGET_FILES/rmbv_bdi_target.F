      module mod_rmbv_bdi
! **********************************************************************
!     Author : C. Voemel
!
!     Date of last modification : 7.7.00
!      
!     Description : PERFORMS MV MULT. WITH MATRIX IN 'BDI'-STORAGE
!                   rmbv = Right Multiplication By Vector: y=Ax
! **********************************************************************

      use representation_of_data
      use properties
      use mod_dense_mat_algos

      implicit none
      
      interface rmbv_bdi
        module procedure irmbv_bdi
        module procedure srmbv_bdi
        module procedure drmbv_bdi
        module procedure crmbv_bdi
        module procedure zrmbv_bdi
      end interface
      
      contains
! **********************************************************************
!
! **********************************************************************
#include "integerrmbv_bdi_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singlePrecisionrmbv_bdi_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doublePrecisionrmbv_bdi_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singleComplexrmbv_bdi_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doubleComplexrmbv_bdi_source.F"
! **********************************************************************
!
! **********************************************************************
      end module mod_rmbv_bdi
