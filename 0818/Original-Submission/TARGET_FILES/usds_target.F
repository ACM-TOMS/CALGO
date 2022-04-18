      module mod_usds
! **********************************************************************
!     Author : C. Voemel
!
!     Date of last modification : 7.7.00
!      
!     Description : RELEASES HANDLES, LOOKS FOR THE "FREEDOM OF MEMORY"
! **********************************************************************

      use representation_of_data
      use properties
      implicit none

      contains

      subroutine usds(nmb,ierr)
      implicit none
      intrinsic modulo
      integer, intent(in) :: nmb
      integer, intent(out) :: ierr
      type(ispmat),pointer :: isp_data
      type(sspmat),pointer :: ssp_data
      type(dspmat),pointer :: dsp_data
      type(cspmat),pointer :: csp_data
      type(zspmat),pointer :: zsp_data
      integer :: rest,val
      
      rest = modulo(nmb,no_of_types)

      select case(rest)
      case(ISP_MATRIX)
! **********************************************************************
#include "integerusds_source.F"
! **********************************************************************
      case(SSP_MATRIX)
! **********************************************************************
#include "singlePrecisionusds_source.F"
! **********************************************************************
      case(DSP_MATRIX)
! **********************************************************************
#include "doublePrecisionusds_source.F"
! **********************************************************************
      case(CSP_MATRIX)
! **********************************************************************
#include "singleComplexusds_source.F"
! **********************************************************************
      case(ZSP_MATRIX)
! **********************************************************************
#include "doubleComplexusds_source.F"
! **********************************************************************
      case default
         ierr = blas_error_param 
      end select
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      end subroutine usds
      end module mod_usds
