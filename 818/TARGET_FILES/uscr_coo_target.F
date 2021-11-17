      module mod_uscr_coo
! **********************************************************************
!     Author : C. Voemel
!
!     Date of last modification : 7.7.00
!      
!     Description : CREATION ROUTINE FOR MATRIX HANDLE FROM 'COO'-FORMAT
! **********************************************************************

      use representation_of_data
      use properties
      use mod_usds

      implicit none
      
      interface uscr_coo
        module procedure iuscr_coo
        module procedure suscr_coo
        module procedure duscr_coo
        module procedure cuscr_coo
        module procedure zuscr_coo
      end interface
      
      contains
! **********************************************************************
!
! **********************************************************************
#include "integeruscr_coo_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singlePrecisionuscr_coo_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doublePrecisionuscr_coo_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singleComplexuscr_coo_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doubleComplexuscr_coo_source.F"
! **********************************************************************
!
! **********************************************************************
      end module mod_uscr_coo
