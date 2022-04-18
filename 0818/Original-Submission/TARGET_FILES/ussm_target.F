      module mod_ussm
! **********************************************************************
!     Author : C. Voemel
!
!     Date of last modification : 7.7.00
!      
!     Description : TRI. SOLVE, CHOOSES APPROPRIATE SUBROUTINES
! **********************************************************************

      use representation_of_data
      use properties
      use mod_sbv

      implicit none

      interface ussm
        module procedure iussm
        module procedure sussm
        module procedure dussm
        module procedure cussm
        module procedure zussm
      end interface

      contains
! **********************************************************************
!
! **********************************************************************
#include "integerussm_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singlePrecisionussm_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doublePrecisionussm_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singleComplexussm_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doubleComplexussm_source.F"
! **********************************************************************
!
! **********************************************************************
      end module  mod_ussm
