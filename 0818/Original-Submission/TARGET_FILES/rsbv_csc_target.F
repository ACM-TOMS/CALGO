      module mod_rsbv_csc
! **********************************************************************
!     Author : C. Voemel
!
!     Date of last modification : 7.7.00
!      
!     Description : PERFORMS TRI. SOLVE WITH MATRIX IN 'CSC'-STORAGE
!                   rsbv = Right Solve By Vector
! **********************************************************************
      use representation_of_data
      use properties

      implicit none

      interface rsbv_csc
        module procedure irsbv_csc
        module procedure srsbv_csc
        module procedure drsbv_csc
        module procedure crsbv_csc
        module procedure zrsbv_csc
      end interface

      contains
! **********************************************************************
!
! **********************************************************************
#include "integerrsbv_csc_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singlePrecisionrsbv_csc_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doublePrecisionrsbv_csc_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singleComplexrsbv_csc_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doubleComplexrsbv_csc_source.F"
! **********************************************************************
!
! **********************************************************************
      end module mod_rsbv_csc
