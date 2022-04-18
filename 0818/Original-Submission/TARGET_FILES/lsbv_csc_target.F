      module mod_lsbv_csc
! **********************************************************************
!     Author : C. Voemel
!
!     Date of last modification : 7.7.00
!      
!     Description : PERFORMS TRI. SOLVE WITH TRANSPOSE IN 'CSC'-STORAGE
!                   lsbv = Left Solve By Vector
! **********************************************************************      
      use representation_of_data
      use properties
      
      implicit none

      interface lsbv_csc
        module procedure ilsbv_csc
        module procedure slsbv_csc
        module procedure dlsbv_csc
        module procedure clsbv_csc
        module procedure zlsbv_csc
      end interface

      contains
! **********************************************************************
!
! **********************************************************************
#include "integerlsbv_csc_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singlePrecisionlsbv_csc_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doublePrecisionlsbv_csc_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singleComplexlsbv_csc_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doubleComplexlsbv_csc_source.F"
! **********************************************************************
!
! **********************************************************************
      end module mod_lsbv_csc
