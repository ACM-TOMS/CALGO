      module mod_usconv_bco2bdi
      use properties
      use mod_conv_tools
      use representation_of_data
      
      contains

      subroutine usconv_bco2bdi(a,ierr) 
      implicit none
      integer,intent(inout) :: a 
      integer ,intent(inout)::ierr
      integer :: res,BLDA,BNDIAG,mb,kb,lb,rest   
      type(ispmat),pointer :: isp_data
      type(sspmat),pointer :: ssp_data
      type(dspmat),pointer :: dsp_data
      type(cspmat),pointer :: csp_data
      type(zspmat),pointer :: zsp_data
      intrinsic min   

      ierr=-1
      rest = modulo(a,no_of_types)

      select case(rest)
      case(ISP_MATRIX)
! **********************************************************************
#include "integerusconv_bco2bdi_source.F"
! **********************************************************************
      case(SSP_MATRIX)
! **********************************************************************
#include "singlePrecisionusconv_bco2bdi_source.F"
! **********************************************************************
      case(DSP_MATRIX)
! **********************************************************************
#include "doublePrecisionusconv_bco2bdi_source.F"
! **********************************************************************
      case(CSP_MATRIX)
! **********************************************************************
#include "singleComplexusconv_bco2bdi_source.F"
! **********************************************************************
      case(ZSP_MATRIX)
! **********************************************************************
#include "doubleComplexusconv_bco2bdi_source.F"
! **********************************************************************
      case default
         ierr = blas_error_param
         return
      end select
      end subroutine usconv_bco2bdi
      end module mod_usconv_bco2bdi
