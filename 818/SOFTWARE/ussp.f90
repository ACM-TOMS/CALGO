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
  call access_matrix(ipmatrix ,a,istat) 
  ipmatrix %property=m
! **********************************************************************
!!***************************************************************************
      case(SSP_MATRIX)
! **********************************************************************
  call access_matrix(spmatrix ,a,istat) 
  spmatrix %property=m
! **********************************************************************
!!***************************************************************************
      case(DSP_MATRIX)
! **********************************************************************
  call access_matrix(dpmatrix ,a,istat) 
  dpmatrix %property=m
! **********************************************************************
!!*************************************************************************** 
      case(CSP_MATRIX)
! **********************************************************************
  call access_matrix(cpmatrix ,a,istat) 
  cpmatrix %property=m
! **********************************************************************
!!*************************************************************************** 
      case(ZSP_MATRIX)
! **********************************************************************
  call access_matrix(zpmatrix ,a,istat) 
  zpmatrix %property=m
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
      case default
         istat = blas_error_param
         return
      end select
      end subroutine ussp
      end module mod_ussp
