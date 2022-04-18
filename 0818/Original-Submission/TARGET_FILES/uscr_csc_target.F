      module mod_uscr_csc
! **********************************************************************
!     Author : C. Voemel
!
!     Date of last modification : 7.7.00
!      
!     Description : CREATION ROUTINE FOR MATRIX HANDLE FROM 'CSC'-FORMAT
! **********************************************************************

      use representation_of_data
      use properties
      use mod_usds
      
      implicit none
      
      interface uscr_csc
        module procedure iuscr_csc
        module procedure suscr_csc
        module procedure duscr_csc
        module procedure cuscr_csc
        module procedure zuscr_csc
      end interface
      
      contains
! **********************************************************************
!
! **********************************************************************
#include "integeruscr_csc_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singlePrecisionuscr_csc_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doublePrecisionuscr_csc_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singleComplexuscr_csc_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doubleComplexuscr_csc_source.F"
! **********************************************************************
!
! **********************************************************************
      end module mod_uscr_csc
