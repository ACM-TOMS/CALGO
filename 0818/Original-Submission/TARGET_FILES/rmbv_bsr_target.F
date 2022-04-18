      module mod_rmbv_bsr
! **********************************************************************
!     Author : C. Voemel
!
!     Date of last modification : 7.7.00
!      
!     Description : PERFORMS MV MULT. WITH MATRIX IN 'BSR'-STORAGE
!                   rmbv = Right Multiplication By Vector: y=Ax
! **********************************************************************

      use representation_of_data
      use properties
      use mod_dense_mat_algos

      implicit none
      
      interface rmbv_bsr
        module procedure irmbv_bsr
        module procedure srmbv_bsr
        module procedure drmbv_bsr
        module procedure crmbv_bsr
        module procedure zrmbv_bsr
      end interface
      
      contains
! **********************************************************************
!
! **********************************************************************
#include "integerrmbv_bsr_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singlePrecisionrmbv_bsr_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doublePrecisionrmbv_bsr_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singleComplexrmbv_bsr_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doubleComplexrmbv_bsr_source.F"
! **********************************************************************
!
! **********************************************************************
      end module mod_rmbv_bsr
