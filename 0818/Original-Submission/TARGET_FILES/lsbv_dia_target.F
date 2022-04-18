      module mod_lsbv_dia
! **********************************************************************
!     Author : C. Voemel
!
!     Date of last modification : 7.7.00
!      
!     Description : PERFORMS TRI. SOLVE WITH MATRIX IN 'DIA'-STORAGE
!                   lsbv = Left Solve By Vector
! **********************************************************************
      use representation_of_data
      use properties

      implicit none

      interface lsbv_dia
        module procedure ilsbv_dia
        module procedure slsbv_dia
        module procedure dlsbv_dia
        module procedure clsbv_dia
        module procedure zlsbv_dia
      end interface

      contains
! **********************************************************************
!
! **********************************************************************
#include "integerlsbv_dia_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singlePrecisionlsbv_dia_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doublePrecisionlsbv_dia_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singleComplexlsbv_dia_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doubleComplexlsbv_dia_source.F"
! **********************************************************************
!
! **********************************************************************
      end module mod_lsbv_dia
