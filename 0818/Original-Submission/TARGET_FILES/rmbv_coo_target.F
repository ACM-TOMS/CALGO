      module mod_rmbv_coo
! **********************************************************************
!     Author : C. Voemel
!
!     Date of last modification : 7.7.00
!      
!     Description : PERFORMS MV MULT. WITH MATRIX IN 'COO'-STORAGE
!                   rmbv = Right Multiplication By Vector: y=Ax
! **********************************************************************

      use representation_of_data
      use properties
      implicit none
      
      interface rmbv_coo
        module procedure irmbv_coo
        module procedure srmbv_coo
        module procedure drmbv_coo
        module procedure crmbv_coo
        module procedure zrmbv_coo
      end interface
      
      contains
! **********************************************************************
!
! **********************************************************************
#include "integerrmbv_coo_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singlePrecisionrmbv_coo_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doublePrecisionrmbv_coo_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singleComplexrmbv_coo_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doubleComplexrmbv_coo_source.F"
! **********************************************************************
!
! **********************************************************************
      end module mod_rmbv_coo
