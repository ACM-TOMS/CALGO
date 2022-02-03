      module mod_rmbv_bco
! **********************************************************************
!     Author : C. Voemel
!
!     Date of last modification : 7.7.00
!      
!     Description : PERFORMS MV MULT. WITH MATRIX IN 'BCO'-STORAGE
!                   rmbv = Right Multiplication By Vector: y=Ax
! **********************************************************************

      use representation_of_data
      use properties
      use mod_dense_mat_algos

      implicit none
      
      interface rmbv_bco
        module procedure irmbv_bco
        module procedure srmbv_bco
        module procedure drmbv_bco
        module procedure crmbv_bco
        module procedure zrmbv_bco
      end interface
      
      contains
! **********************************************************************
!
! **********************************************************************
#include "integerrmbv_bco_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singlePrecisionrmbv_bco_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doublePrecisionrmbv_bco_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singleComplexrmbv_bco_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doubleComplexrmbv_bco_source.F"
! **********************************************************************
!
! **********************************************************************
      end module mod_rmbv_bco
