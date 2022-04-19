      module mod_uscr_bsc
! **********************************************************************
!     Author : C. Voemel
!
!     Date of last modification : 7.7.00
!      
!     Description : CREATION ROUTINE FOR MATRIX HANDLE FROM 'BSC'-FORMAT
! **********************************************************************

      use representation_of_data
      use properties
      use mod_usds

      implicit none
      
      interface uscr_bsc
        module procedure iuscr_bsc
        module procedure suscr_bsc
        module procedure duscr_bsc
        module procedure cuscr_bsc
        module procedure zuscr_bsc
      end interface
      
      contains
! **********************************************************************
!
! **********************************************************************
#include "integeruscr_bsc_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singlePrecisionuscr_bsc_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doublePrecisionuscr_bsc_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singleComplexuscr_bsc_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doubleComplexuscr_bsc_source.F"
! **********************************************************************
!
! **********************************************************************
      end module mod_uscr_bsc
