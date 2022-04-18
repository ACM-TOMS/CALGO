module mod_uscr_vbr
! **********************************************************************
!     Author : C. Voemel
!
!     Date of last modification : 7.7.00
!      
!     Description : CREATION ROUTINE FOR MATRIX HANDLE FROM 'VBR'-FORMAT
! **********************************************************************

use representation_of_data
use properties
use mod_usds

implicit none

interface uscr_vbr
  module procedure iuscr_vbr
  module procedure suscr_vbr
  module procedure duscr_vbr
  module procedure cuscr_vbr
  module procedure zuscr_vbr
end interface

contains
! **********************************************************************
!
! **********************************************************************
#include "integeruscr_vbr_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singlePrecisionuscr_vbr_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doublePrecisionuscr_vbr_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singleComplexuscr_vbr_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doubleComplexuscr_vbr_source.F"
! **********************************************************************
!
! **********************************************************************
end module mod_uscr_vbr
