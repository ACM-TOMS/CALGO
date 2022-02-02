      module mod_usgp
      use mod_INSERTING
      use properties
      use representation_of_data
      
      contains
      
      subroutine usgp (a,pname,m)
      implicit none
      integer ,intent(in)::a
      integer,intent(out)::m
      integer,intent(in)::pname
      type(ispmat),pointer :: isp_data
      type(sspmat),pointer :: ssp_data
      type(dspmat),pointer :: dsp_data
      type(cspmat),pointer :: csp_data
      type(zspmat),pointer :: zsp_data
      type(i_matrix),pointer :: imatrix
      type(s_matrix),pointer :: smatrix
      type(d_matrix),pointer :: dmatrix
      type(c_matrix),pointer :: cmatrix
      type(z_matrix),pointer :: zmatrix
      integer ::rest,ierr
      character ::test

      rest = modulo(a,no_of_types)
      select case(rest)
! **********************************************************************
!
! **********************************************************************
      case(ISP_MATRIX)
#include "integerusgp_source.F"
! **********************************************************************
!
! **********************************************************************
      case(SSP_MATRIX)
#include "singlePrecisionusgp_source.F"
! **********************************************************************
!
! **********************************************************************
      case(DSP_MATRIX)
#include "doublePrecisionusgp_source.F"
! **********************************************************************
!
! **********************************************************************
      case(CSP_MATRIX)
#include "singleComplexusgp_source.F"
! **********************************************************************
!
! **********************************************************************
      case(ZSP_MATRIX)
#include "doubleComplexusgp_source.F"
! **********************************************************************
!
! **********************************************************************
      case default
         return
      end select

      end subroutine usgp

      end module mod_usgp
