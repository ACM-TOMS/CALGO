      module mod_uscr_bsr
! **********************************************************************
!     Author : C. Voemel
!
!     Date of last modification : 7.7.00
!      
!     Description : CREATION ROUTINE FOR MATRIX HANDLE FROM 'BSR'-FORMAT
! **********************************************************************

      use representation_of_data
      use properties
      use mod_usds

      implicit none
      
      interface uscr_bsr
        module procedure iuscr_bsr
        module procedure suscr_bsr
        module procedure duscr_bsr
        module procedure cuscr_bsr
        module procedure zuscr_bsr
      end interface
      
      contains
! **********************************************************************
!
! **********************************************************************
#include "integeruscr_bsr_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singlePrecisionuscr_bsr_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doublePrecisionuscr_bsr_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singleComplexuscr_bsr_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doubleComplexuscr_bsr_source.F"
! **********************************************************************
!
! **********************************************************************
      end module mod_uscr_bsr
