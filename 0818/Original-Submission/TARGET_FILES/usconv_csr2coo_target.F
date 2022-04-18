      module mod_usconv_csr2coo
      use properties
      use mod_conv_tools
      use representation_of_data
      
      contains
      subroutine usconv_csr2coo(a,ierr) 
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
#include "integerusconv_csr2coo_source.F"
! **********************************************************************
      case(SSP_MATRIX)
! **********************************************************************
#include "singlePrecisionusconv_csr2coo_source.F"
! **********************************************************************
      case(DSP_MATRIX)
! **********************************************************************
#include "doublePrecisionusconv_csr2coo_source.F"
! **********************************************************************
      case(CSP_MATRIX)
! **********************************************************************
#include "singleComplexusconv_csr2coo_source.F"
! **********************************************************************
      case(ZSP_MATRIX)
! **********************************************************************
#include "doubleComplexusconv_csr2coo_source.F"
! **********************************************************************
      case default
         ierr = blas_error_param
         return
      end select
      end subroutine usconv_csr2coo
      end module mod_usconv_csr2coo
