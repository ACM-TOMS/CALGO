      module mod_ussp
      use mod_INSERTING
      use properties
      
      contains
      subroutine ussp(a,m,istat)
      implicit none
      integer ,intent(inout)::a
      integer,intent(in)::m
      integer, intent(out)::istat
      integer::b,rest
      type(i_matrix),pointer ::ipmatrix
      type(d_matrix),pointer ::dpmatrix
      type(s_matrix),pointer ::spmatrix
      type(c_matrix),pointer ::cpmatrix
      type(z_matrix),pointer ::zpmatrix
      
      b=-a
      istat = 0
      rest = modulo(b,no_of_types)
      select case(rest)
      case(ISP_MATRIX)
! **********************************************************************
#include "integerussp_source.F"
! **********************************************************************
!!***************************************************************************
      case(SSP_MATRIX)
! **********************************************************************
#include "singlePrecisionussp_source.F"
! **********************************************************************
!!***************************************************************************
      case(DSP_MATRIX)
! **********************************************************************
#include "doublePrecisionussp_source.F"
! **********************************************************************
!!*************************************************************************** 
      case(CSP_MATRIX)
! **********************************************************************
#include "singleComplexussp_source.F"
! **********************************************************************
!!*************************************************************************** 
      case(ZSP_MATRIX)
! **********************************************************************
#include "doubleComplexussp_source.F"
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
      case default
         istat = blas_error_param
         return
      end select
      end subroutine ussp
      end module mod_ussp
