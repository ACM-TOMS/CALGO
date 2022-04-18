module mod_uscr_bco
! **********************************************************************
!     Author : C. Voemel
!
!     Date of last modification : 7.7.00
!      
!     Description : CREATION ROUTINE FOR MATRIX HANDLE FROM 'BCO'-FORMAT
! **********************************************************************

use representation_of_data
use properties
use mod_usds

implicit none

interface uscr_bco
  module procedure iuscr_bco
  module procedure suscr_bco
  module procedure duscr_bco
  module procedure cuscr_bco
  module procedure zuscr_bco
end interface

contains
! **********************************************************************
!
! **********************************************************************
#include "integeruscr_bco_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singlePrecisionuscr_bco_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doublePrecisionuscr_bco_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singleComplexuscr_bco_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doubleComplexuscr_bco_source.F"
! **********************************************************************
!
! **********************************************************************
end module mod_uscr_bco
