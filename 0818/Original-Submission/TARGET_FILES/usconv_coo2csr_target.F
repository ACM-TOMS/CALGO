      module mod_usconv_coo2csr
      use properties
      use mod_conv_tools
      use representation_of_data
      
      
      contains
      subroutine usconv_coo2csr(a,ierr) 
      implicit none
      integer,intent(inout) :: a 
      type(ispmat),pointer :: isp_data
      type(sspmat),pointer :: ssp_data
      type(dspmat),pointer :: dsp_data
      type(cspmat),pointer :: csp_data
      type(zspmat),pointer :: zsp_data
      integer ,intent(inout)::ierr
      integer :: res,rest
     
      ierr=-1
      rest = modulo(a,no_of_types)
      select case(rest)
      case(ISP_MATRIX) 
!!*************************************************************************** 
! **********************************************************************
#include "integerusconv_coo2csr_source.F"
      case(SSP_MATRIX) 
!!*************************************************************************** 
! **********************************************************************
#include "singlePrecisionusconv_coo2csr_source.F"
! **********************************************************************
!!***************************************************************************
      case(DSP_MATRIX) 
! **********************************************************************
#include "doublePrecisionusconv_coo2csr_source.F"
! **********************************************************************
!!***************************************************************************
      case(CSP_MATRIX) 
! **********************************************************************
#include "singleComplexusconv_coo2csr_source.F"
! **********************************************************************
!!***************************************************************************
      case(ZSP_MATRIX) 
! **********************************************************************
#include "doubleComplexusconv_coo2csr_source.F"
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************  
      case default
         ierr = blas_error_param
         return
      end select
      end subroutine usconv_coo2csr
      end module mod_usconv_coo2csr
