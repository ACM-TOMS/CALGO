      module mod_ussv
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

      interface ussv
        module procedure iussv
        module procedure sussv
        module procedure dussv
        module procedure cussv
        module procedure zussv
      end interface

      contains
! **********************************************************************
!
! **********************************************************************
#include "integerussv_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singlePrecisionussv_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doublePrecisionussv_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singleComplexussv_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doubleComplexussv_source.F"
! **********************************************************************
!
! **********************************************************************
      end module  mod_ussv
