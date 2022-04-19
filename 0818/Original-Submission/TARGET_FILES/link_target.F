      module  representation_of_data
! **********************************************************************
!     Author : C. Voemel
!
!     Date of last modification : 7.7.00
!      
!     Description : THE PRINCIPAL DATA STRUCTURE
!                   Matrix data is stored in nodes of a linked list
!                   Node number is the handle number
!                   new: creates new node WITHOUT initialization
!                   del: frees unused memory, does NOT care if there
!                        is other memory that should be freed first
!                   accessdata: given a handle number, it returns a
!                        pointer to the matrix inside the relevant node
!
! **********************************************************************

      use types
      use properties

      implicit none
      
      interface accessdata 
         module procedure accessdata_isp
         module procedure accessdata_ssp
         module procedure accessdata_dsp
         module procedure accessdata_csp
         module procedure accessdata_zsp
      end interface

      type isp_linknode
        type(ispmat) :: contents
        integer :: number
        type(isp_linknode), pointer :: pntr
      end type isp_linknode
      type ssp_linknode
        type(sspmat) :: contents
        integer :: number
        type(ssp_linknode), pointer :: pntr
      end type ssp_linknode
      type dsp_linknode
        type(dspmat) :: contents
        integer :: number
        type(dsp_linknode), pointer :: pntr
      end type dsp_linknode
      type csp_linknode
        type(cspmat) :: contents
        integer :: number
        type(csp_linknode), pointer :: pntr
      end type csp_linknode
      type zsp_linknode
        type(zspmat) :: contents
        integer :: number
        type(zsp_linknode), pointer :: pntr
      end type zsp_linknode
 
      type(isp_linknode), pointer,SAVE,PRIVATE :: isp_first, isp_last
      type(ssp_linknode), pointer,SAVE,PRIVATE :: ssp_first, ssp_last
      type(dsp_linknode), pointer,SAVE,PRIVATE :: dsp_first, dsp_last
      type(csp_linknode), pointer,SAVE,PRIVATE :: csp_first, csp_last
      type(zsp_linknode), pointer,SAVE,PRIVATE :: zsp_first, zsp_last

      logical,SAVE,PRIVATE :: isp_init = .FALSE.
      logical,SAVE,PRIVATE :: ssp_init = .FALSE.
      logical,SAVE,PRIVATE :: dsp_init = .FALSE.
      logical,SAVE,PRIVATE :: csp_init = .FALSE.
      logical,SAVE,PRIVATE :: zsp_init = .FALSE.

      contains 
! **********************************************************************
!
! **********************************************************************
#include "integerlink_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singlePrecisionlink_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doublePrecisionlink_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singleComplexlink_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doubleComplexlink_source.F"
! **********************************************************************
!
! **********************************************************************
      end module representation_of_data
