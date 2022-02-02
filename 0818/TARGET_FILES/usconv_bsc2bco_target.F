      module mod_usconv_bsc2bco
      use properties
      use mod_conv_tools
      use representation_of_data

      contains
      subroutine usconv_bsc2bco(a,ierr) 
      implicit none
      integer,intent(inout) :: a 
      integer ,intent(inout)::ierr
      integer :: res,s,rest
      type(ispmat),pointer :: isp_data
      type(sspmat),pointer :: ssp_data
      type(dspmat),pointer :: dsp_data
      type(cspmat),pointer :: csp_data
      type(zspmat),pointer :: zsp_data
      
      ierr=-1
      rest = modulo(a,no_of_types)
      
      select case(rest)
      case(ISP_MATRIX)
! **********************************************************************
#include "integerusconv_bsc2bco_source.F"
! **********************************************************************
      case(SSP_MATRIX)
! **********************************************************************
#include "singlePrecisionusconv_bsc2bco_source.F"
! **********************************************************************
      case(DSP_MATRIX)
! **********************************************************************
#include "doublePrecisionusconv_bsc2bco_source.F"
! **********************************************************************
      case(CSP_MATRIX)
! **********************************************************************
#include "singleComplexusconv_bsc2bco_source.F"
! **********************************************************************
      case(ZSP_MATRIX)
! **********************************************************************
#include "doubleComplexusconv_bsc2bco_source.F"
! **********************************************************************
      case default
         ierr = blas_error_param
         return
      end select
      end subroutine usconv_bsc2bco
      end module mod_usconv_bsc2bco
