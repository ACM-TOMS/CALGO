      module mod_uscr_dia
! **********************************************************************
!     Author : C. Voemel
!
!     Date of last modification : 7.7.00
!      
!     Description : CREATION ROUTINE FOR MATRIX HANDLE FROM 'DIA'-FORMAT
! **********************************************************************

      use representation_of_data
      use properties

      implicit none
      
      interface uscr_dia
        module procedure iuscr_dia
        module procedure suscr_dia
        module procedure duscr_dia
        module procedure cuscr_dia
        module procedure zuscr_dia
      end interface
      
      contains
! **********************************************************************
!
! **********************************************************************
#include "integeruscr_dia_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singlePrecisionuscr_dia_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doublePrecisionuscr_dia_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singleComplexuscr_dia_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doubleComplexuscr_dia_source.F"
! **********************************************************************
!
! **********************************************************************
      end module mod_uscr_dia
