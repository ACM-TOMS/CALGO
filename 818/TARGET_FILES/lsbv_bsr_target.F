      module mod_lsbv_bsr
! **********************************************************************
!     Author : C. Voemel
!
!     Date of last modification : 7.7.00
!      
!     Description : PERFORMS TRI. SOLVE WITH TRANSPOSE IN 'BSR'-STORAGE
!                   lsbv = Left Solve By Vector
! **********************************************************************      
      use representation_of_data
      use properties
      use mod_dense_mat_algos
      
      implicit none

      interface lsbv_bsr
        module procedure ilsbv_bsr
        module procedure slsbv_bsr
        module procedure dlsbv_bsr
        module procedure clsbv_bsr
        module procedure zlsbv_bsr
      end interface

      contains
! **********************************************************************
!
! **********************************************************************
#include "integerlsbv_bsr_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singlePrecisionlsbv_bsr_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doublePrecisionlsbv_bsr_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singleComplexlsbv_bsr_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doubleComplexlsbv_bsr_source.F"
! **********************************************************************
!
! **********************************************************************
      end module mod_lsbv_bsr
