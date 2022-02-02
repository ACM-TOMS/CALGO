      module mod_info
! **********************************************************************
!     Author : C. Voemel
!
!     Date of last modification : 7.7.00
!      
!     Description : FOR DEBUGGING ONLY !!! 
!                   "print" displays data for given handle number
! **********************************************************************

      use representation_of_data
      use properties

      implicit none

      contains

      subroutine print(nmb,ierr)
      implicit none
      intrinsic modulo
      integer, intent(in) :: nmb
      integer, intent(out) :: ierr
      type(ispmat),pointer :: isp_data
      type(sspmat),pointer :: ssp_data
      type(dspmat),pointer :: dsp_data
      type(cspmat),pointer :: csp_data
      type(zspmat),pointer :: zsp_data
      integer :: rest,base,copy,nnz,rowdim,coldim
      character :: style,diag,type,part,store

      rest = modulo(nmb,no_of_types)

      select case(rest)
      case(ISP_MATRIX)
! **********************************************************************
#include "integerinfo_source.F"
! **********************************************************************
      case(SSP_MATRIX)
! **********************************************************************
#include "singlePrecisioninfo_source.F"
! **********************************************************************
      case(DSP_MATRIX)
! **********************************************************************
#include "doublePrecisioninfo_source.F"
! **********************************************************************
      case(CSP_MATRIX)
! **********************************************************************
#include "singleComplexinfo_source.F"
! **********************************************************************
      case(ZSP_MATRIX)
! **********************************************************************
#include "doubleComplexinfo_source.F"
! **********************************************************************
      case default
         write(*,*) 'Wrong matrix type !'
         ierr = -1 
      end select
      
      end subroutine print
      end module mod_info
