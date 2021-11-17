      module mod_rmbv_bsc
! **********************************************************************
!     Author : C. Voemel
!
!     Date of last modification : 7.7.00
!      
!     Description : PERFORMS MV MULT. WITH MATRIX IN 'BSC'-STORAGE
!                   rmbv = Right Multiplication By Vector: y=Ax
! **********************************************************************

      use representation_of_data
      use properties
      use mod_dense_mat_algos

      implicit none
      
      interface rmbv_bsc
        module procedure irmbv_bsc
        module procedure srmbv_bsc
        module procedure drmbv_bsc
        module procedure crmbv_bsc
        module procedure zrmbv_bsc
      end interface
      
      contains
! **********************************************************************
!
! **********************************************************************
#include "integerrmbv_bsc_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singlePrecisionrmbv_bsc_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doublePrecisionrmbv_bsc_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singleComplexrmbv_bsc_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doubleComplexrmbv_bsc_source.F"
! **********************************************************************
!
! **********************************************************************
      end module mod_rmbv_bsc
