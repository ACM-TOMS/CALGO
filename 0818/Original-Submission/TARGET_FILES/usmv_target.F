      module mod_usmv
! **********************************************************************
!     Author : C. Voemel
!
!     Date of last modification : 7.7.00
!      
!     Description : MV MULTIPLICATION, CHOOSES APPROPRIATE SUBROUTINES
! **********************************************************************

      use representation_of_data
      use properties
      use mod_mbv

      implicit none

      interface usmv
        module procedure iusmv
        module procedure susmv
        module procedure dusmv
        module procedure cusmv
        module procedure zusmv
      end interface

      contains
! **********************************************************************
!
! **********************************************************************
#include "integerusmv_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singlePrecisionusmv_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doublePrecisionusmv_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singleComplexusmv_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doubleComplexusmv_source.F"
! **********************************************************************
!
! **********************************************************************
      end module mod_usmv
