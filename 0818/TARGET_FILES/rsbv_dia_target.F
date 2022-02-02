      module mod_rsbv_dia
! **********************************************************************
!     Author : C. Voemel
!
!     Date of last modification : 7.7.00
!      
!     Description : PERFORMS TRI. SOLVE WITH MATRIX IN 'DIA'-STORAGE
!                   rsbv = Right Solve By Vector
! **********************************************************************
      use representation_of_data
      use properties

      implicit none

      interface rsbv_dia
        module procedure irsbv_dia
        module procedure srsbv_dia
        module procedure drsbv_dia
        module procedure crsbv_dia
        module procedure zrsbv_dia
      end interface

      contains
! **********************************************************************
!
! **********************************************************************
#include "integerrsbv_dia_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singlePrecisionrsbv_dia_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doublePrecisionrsbv_dia_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singleComplexrsbv_dia_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doubleComplexrsbv_dia_source.F"
! **********************************************************************
!
! **********************************************************************
      end module mod_rsbv_dia
