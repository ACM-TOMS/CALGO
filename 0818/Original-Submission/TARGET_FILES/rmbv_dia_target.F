      module mod_rmbv_dia
! **********************************************************************
!     Author : C. Voemel
!
!     Date of last modification : 7.7.00
!      
!     Description : PERFORMS MV MULT. WITH MATRIX IN 'DIA'-STORAGE
!                   rmbv = Right Multiplication By Vector: y=Ax
! **********************************************************************

      use representation_of_data
      use properties
      implicit none
      
      interface rmbv_dia
        module procedure irmbv_dia
        module procedure srmbv_dia
        module procedure drmbv_dia
        module procedure crmbv_dia
        module procedure zrmbv_dia
      end interface
      
      contains
! **********************************************************************
!
! **********************************************************************
#include "integerrmbv_dia_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singlePrecisionrmbv_dia_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doublePrecisionrmbv_dia_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singleComplexrmbv_dia_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doubleComplexrmbv_dia_source.F"
! **********************************************************************
!
! **********************************************************************
      end module mod_rmbv_dia
