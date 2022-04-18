      module mod_usmm
! **********************************************************************
!     Author : C. Voemel
!
!     Date of last modification : 7.7.00
!      
!     Description : MM MULTIPLICATION, CHOOSES APPROPRIATE SUBROUTINES
! **********************************************************************

      use representation_of_data
      use properties
      use mod_mbv

      implicit none

      interface usmm
        module procedure iusmm
        module procedure susmm
        module procedure dusmm
        module procedure cusmm
        module procedure zusmm
      end interface

      contains
! **********************************************************************
!
! **********************************************************************
#include "integerusmm_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singlePrecisionusmm_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doublePrecisionusmm_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singleComplexusmm_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doubleComplexusmm_source.F"
! **********************************************************************
!
! **********************************************************************
      end module mod_usmm
