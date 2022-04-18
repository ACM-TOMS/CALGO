module mod_uscr_insert_col
use blas_sparse_namedconstants
use mod_uscr_insert_entry
interface uscr_insert_col
   module procedure iuscr_insert_col
   module procedure suscr_insert_col
   module procedure duscr_insert_col
   module procedure cuscr_insert_col
   module procedure zuscr_insert_col
end interface
contains
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
      subroutine  iuscr_insert_col  (a,j,val,indx,istat)
      implicit none
      integer  ,dimension(:),intent(in) ::val
      integer ,intent(in) ::a,j
      integer ,intent(out) ::istat
      integer,dimension(:),intent(in)::indx
      integer ::i,s
      istat=-1
      s=size(val)
      do i=1,s
         call   iuscr_insert_entry (a,val(i),indx(i),j,istat)
         if(istat.ne.0) return
      end do
      end subroutine  iuscr_insert_col  
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
      subroutine  suscr_insert_col  (a,j,val,indx,istat)
      implicit none
      real(KIND=sp)  ,dimension(:),intent(in) ::val
      integer ,intent(in) ::a,j
      integer ,intent(out) ::istat
      integer,dimension(:),intent(in)::indx
      integer ::i,s
      istat=-1
      s=size(val)
      do i=1,s
         call   suscr_insert_entry (a,val(i),indx(i),j,istat)
         if(istat.ne.0) return
      end do
      end subroutine  suscr_insert_col  
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
      subroutine  duscr_insert_col  (a,j,val,indx,istat)
      implicit none
      real(KIND=dp)  ,dimension(:),intent(in) ::val
      integer ,intent(in) ::a,j
      integer ,intent(out) ::istat
      integer,dimension(:),intent(in)::indx
      integer ::i,s
      istat=-1
      s=size(val)
      do i=1,s
         call   duscr_insert_entry (a,val(i),indx(i),j,istat)
         if(istat.ne.0) return
      end do
      end subroutine  duscr_insert_col  
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
      subroutine  cuscr_insert_col  (a,j,val,indx,istat)
      implicit none
      complex(KIND=sp)  ,dimension(:),intent(in) ::val
      integer ,intent(in) ::a,j
      integer ,intent(out) ::istat
      integer,dimension(:),intent(in)::indx
      integer ::i,s
      istat=-1
      s=size(val)
      do i=1,s
         call   cuscr_insert_entry (a,val(i),indx(i),j,istat)
         if(istat.ne.0) return
      end do
      end subroutine  cuscr_insert_col  
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
      subroutine  zuscr_insert_col  (a,j,val,indx,istat)
      implicit none
      complex(KIND=dp)  ,dimension(:),intent(in) ::val
      integer ,intent(in) ::a,j
      integer ,intent(out) ::istat
      integer,dimension(:),intent(in)::indx
      integer ::i,s
      istat=-1
      s=size(val)
      do i=1,s
         call   zuscr_insert_entry (a,val(i),indx(i),j,istat)
         if(istat.ne.0) return
      end do
      end subroutine  zuscr_insert_col  
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
end module mod_uscr_insert_col
