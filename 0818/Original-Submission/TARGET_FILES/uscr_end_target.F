      module mod_uscr_end
      use mod_INSERTING
      use mod_INS_ROUTINER
      use properties
      
   
      contains

      subroutine uscr_end(a,istat)
      implicit none
      integer ,intent(inout)::a,istat
      integer::prpty,rest,b
      type(i_matrix),pointer ::ipmatrix
      type(d_matrix),pointer ::dpmatrix
      type(s_matrix),pointer ::spmatrix
      type(c_matrix),pointer ::cpmatrix
      type(z_matrix),pointer ::zpmatrix
      
      b=-a
      rest = modulo(b,no_of_types)
      select case(rest)
      case(ISP_MATRIX)
! **********************************************************************
#include "integeruscr_end_source.F"
!!*************************************************************************** 
! **********************************************************************
      case(SSP_MATRIX)
! **********************************************************************
#include "singlePrecisionuscr_end_source.F"
! **********************************************************************
!!***************************************************************************
      case(DSP_MATRIX)
! **********************************************************************
#include "doublePrecisionuscr_end_source.F"
! **********************************************************************
!!***************************************************************************
      case(CSP_MATRIX)
! **********************************************************************
#include "singleComplexuscr_end_source.F"
! **********************************************************************
!!*************************************************************************** 
      case(ZSP_MATRIX)
! **********************************************************************
#include "doubleComplexuscr_end_source.F"
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
      case default
         istat = blas_error_param
         return
      end select
      end subroutine uscr_end
      end module mod_uscr_end
