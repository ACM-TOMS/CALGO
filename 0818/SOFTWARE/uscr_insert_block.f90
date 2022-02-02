module mod_uscr_insert_block
use blas_sparse_namedconstants
use mod_INSERTING
use mod_INS_ROUTINER
interface uscr_insert_block
   module procedure iuscr_insert_block
   module procedure suscr_insert_block
   module procedure duscr_insert_block
   module procedure cuscr_insert_block
   module procedure zuscr_insert_block
end interface
 contains
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
      subroutine  iuscr_insert_block (a,val,bi,bj,istat)
      implicit none
      integer  ,dimension(:,:),intent(in) ::val
      integer ,intent(in) ::a,bi,bj
      integer,intent(out)::istat
      type(i_matrix ),pointer ::pmatrix
      istat=-1
      call access_matrix(pmatrix,a,istat) 
      if(istat.ne.0) return
      select case(pmatrix%format)
      case('block')
         call iINS_block (pmatrix,val,bi,bj,istat)
      case('vblock')
         call  iINS_varblock (pmatrix,val,bi,bj,istat)
      case default
         istat = blas_error_param    
         return
      end select
      end subroutine  iuscr_insert_block 
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
      subroutine  suscr_insert_block (a,val,bi,bj,istat)
      implicit none
      real(KIND=sp)  ,dimension(:,:),intent(in) ::val
      integer ,intent(in) ::a,bi,bj
      integer,intent(out)::istat
      type(s_matrix ),pointer ::pmatrix
      istat=-1
      call access_matrix(pmatrix,a,istat) 
      if(istat.ne.0) return
      select case(pmatrix%format)
      case('block')
         call sINS_block (pmatrix,val,bi,bj,istat)
      case('vblock')
         call  sINS_varblock (pmatrix,val,bi,bj,istat)
      case default
         istat = blas_error_param    
         return
      end select
      end subroutine  suscr_insert_block 
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
      subroutine  duscr_insert_block (a,val,bi,bj,istat)
      implicit none
      real(KIND=dp)  ,dimension(:,:),intent(in) ::val
      integer ,intent(in) ::a,bi,bj
      integer,intent(out)::istat
      type(d_matrix ),pointer ::pmatrix
      istat=-1
      call access_matrix(pmatrix,a,istat) 
      if(istat.ne.0) return
      select case(pmatrix%format)
      case('block')
         call dINS_block (pmatrix,val,bi,bj,istat)
      case('vblock')
         call  dINS_varblock (pmatrix,val,bi,bj,istat)
      case default
         istat = blas_error_param    
         return
      end select
      end subroutine  duscr_insert_block 
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
      subroutine  cuscr_insert_block (a,val,bi,bj,istat)
      implicit none
      complex(KIND=sp)  ,dimension(:,:),intent(in) ::val
      integer ,intent(in) ::a,bi,bj
      integer,intent(out)::istat
      type(c_matrix ),pointer ::pmatrix
      istat=-1
      call access_matrix(pmatrix,a,istat) 
      if(istat.ne.0) return
      select case(pmatrix%format)
      case('block')
         call cINS_block (pmatrix,val,bi,bj,istat)
      case('vblock')
         call  cINS_varblock (pmatrix,val,bi,bj,istat)
      case default
         istat = blas_error_param    
         return
      end select
      end subroutine  cuscr_insert_block 
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
      subroutine  zuscr_insert_block (a,val,bi,bj,istat)
      implicit none
      complex(KIND=dp)  ,dimension(:,:),intent(in) ::val
      integer ,intent(in) ::a,bi,bj
      integer,intent(out)::istat
      type(z_matrix ),pointer ::pmatrix
      istat=-1
      call access_matrix(pmatrix,a,istat) 
      if(istat.ne.0) return
      select case(pmatrix%format)
      case('block')
         call zINS_block (pmatrix,val,bi,bj,istat)
      case('vblock')
         call  zINS_varblock (pmatrix,val,bi,bj,istat)
      case default
         istat = blas_error_param    
         return
      end select
      end subroutine  zuscr_insert_block 
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
end module mod_uscr_insert_block
