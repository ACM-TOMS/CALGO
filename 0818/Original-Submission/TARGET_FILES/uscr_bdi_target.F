      module mod_uscr_bdi
! **********************************************************************
!     Author : C. Voemel
!
!     Date of last modification : 7.7.00
!      
!     Description : CREATION ROUTINE FOR MATRIX HANDLE FROM 'BDI'-FORMAT
! **********************************************************************

      use representation_of_data
      use properties
      use mod_usds

      implicit none
      
      interface uscr_bdi
        module procedure iuscr_bdi
        module procedure suscr_bdi
        module procedure duscr_bdi
        module procedure cuscr_bdi
        module procedure zuscr_bdi
      end interface
      
      contains
! **********************************************************************
!
! **********************************************************************
#include "integeruscr_bdi_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singlePrecisionuscr_bdi_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doublePrecisionuscr_bdi_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singleComplexuscr_bdi_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doubleComplexuscr_bdi_source.F"
! **********************************************************************
!
! **********************************************************************
      end module mod_uscr_bdi
