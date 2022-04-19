module mod_uscr_insert_entry
use blas_sparse_namedconstants
use mod_INSERTING
use mod_INS_ROUTINER
implicit none
interface uscr_insert_entry
module procedure iuscr_insert_entry
module procedure suscr_insert_entry
module procedure duscr_insert_entry
module procedure cuscr_insert_entry
module procedure zuscr_insert_entry
end interface
contains
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
      subroutine iuscr_insert_entry  (a,val,i,j,istat)
      implicit none
      integer  ,intent(in) ::val
      integer ,intent(in) ::a,i,j
      integer,intent(out)::istat
      type(i_matrix ),pointer ::pmatrix
      istat=-1
      call access_matrix(pmatrix,a,istat) 
      if (istat.ne.0) return
      select case(pmatrix%format)
      case('block')
         call iINS_bl_entr (pmatrix,val,i,j,istat)
      case('vblock')
         call  iINS_varbl_entr (pmatrix,val,i,j,istat) 
      case('normal')
         call iINS_entry (pmatrix,val,i,j,istat)
      case default
         istat = blas_error_param
         return
      end select
      end subroutine iuscr_insert_entry  
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
      subroutine suscr_insert_entry  (a,val,i,j,istat)
      implicit none
      real(KIND=sp)  ,intent(in) ::val
      integer ,intent(in) ::a,i,j
      integer,intent(out)::istat
      type(s_matrix ),pointer ::pmatrix
      istat=-1
      call access_matrix(pmatrix,a,istat) 
      if (istat.ne.0) return
      select case(pmatrix%format)
      case('block')
         call sINS_bl_entr (pmatrix,val,i,j,istat)
      case('vblock')
         call  sINS_varbl_entr (pmatrix,val,i,j,istat) 
      case('normal')
         call sINS_entry (pmatrix,val,i,j,istat)
      case default
         istat = blas_error_param
         return
      end select
      end subroutine suscr_insert_entry  
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
      subroutine duscr_insert_entry  (a,val,i,j,istat)
      implicit none
      real(KIND=dp)  ,intent(in) ::val
      integer ,intent(in) ::a,i,j
      integer,intent(out)::istat
      type(d_matrix ),pointer ::pmatrix
      istat=-1
      call access_matrix(pmatrix,a,istat) 
      if (istat.ne.0) return
      select case(pmatrix%format)
      case('block')
         call dINS_bl_entr (pmatrix,val,i,j,istat)
      case('vblock')
         call  dINS_varbl_entr (pmatrix,val,i,j,istat) 
      case('normal')
         call dINS_entry (pmatrix,val,i,j,istat)
      case default
         istat = blas_error_param
         return
      end select
      end subroutine duscr_insert_entry  
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
      subroutine cuscr_insert_entry  (a,val,i,j,istat)
      implicit none
      complex(KIND=sp)  ,intent(in) ::val
      integer ,intent(in) ::a,i,j
      integer,intent(out)::istat
      type(c_matrix ),pointer ::pmatrix
      istat=-1
      call access_matrix(pmatrix,a,istat) 
      if (istat.ne.0) return
      select case(pmatrix%format)
      case('block')
         call cINS_bl_entr (pmatrix,val,i,j,istat)
      case('vblock')
         call  cINS_varbl_entr (pmatrix,val,i,j,istat) 
      case('normal')
         call cINS_entry (pmatrix,val,i,j,istat)
      case default
         istat = blas_error_param
         return
      end select
      end subroutine cuscr_insert_entry  
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
      subroutine zuscr_insert_entry  (a,val,i,j,istat)
      implicit none
      complex(KIND=dp)  ,intent(in) ::val
      integer ,intent(in) ::a,i,j
      integer,intent(out)::istat
      type(z_matrix ),pointer ::pmatrix
      istat=-1
      call access_matrix(pmatrix,a,istat) 
      if (istat.ne.0) return
      select case(pmatrix%format)
      case('block')
         call zINS_bl_entr (pmatrix,val,i,j,istat)
      case('vblock')
         call  zINS_varbl_entr (pmatrix,val,i,j,istat) 
      case('normal')
         call zINS_entry (pmatrix,val,i,j,istat)
      case default
         istat = blas_error_param
         return
      end select
      end subroutine zuscr_insert_entry  
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
end module mod_uscr_insert_entry
