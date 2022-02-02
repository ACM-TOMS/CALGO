      module mod_rmbv_csc
! **********************************************************************
!     Author : C. Voemel
!
!     Date of last modification : 7.7.00
!      
!     Description : PERFORMS MV MULT. WITH MATRIX IN 'CSC'-STORAGE
!                   rmbv = Right Multiplication By Vector: y=Ax
! **********************************************************************

      use representation_of_data
      use properties
      implicit none
      
      interface rmbv_csc
        module procedure irmbv_csc
        module procedure srmbv_csc
        module procedure drmbv_csc
        module procedure crmbv_csc
        module procedure zrmbv_csc
      end interface
      
      contains
! **********************************************************************
!
! **********************************************************************
#include "integerrmbv_csc_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singlePrecisionrmbv_csc_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doublePrecisionrmbv_csc_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singleComplexrmbv_csc_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doubleComplexrmbv_csc_source.F"
! **********************************************************************
!
! **********************************************************************
      end module mod_rmbv_csc
