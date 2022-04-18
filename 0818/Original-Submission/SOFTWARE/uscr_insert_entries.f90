module mod_uscr_insert_entries
use blas_sparse_namedconstants
use mod_uscr_insert_entry
interface uscr_insert_entries
   module procedure iuscr_insert_entries
   module procedure suscr_insert_entries
   module procedure duscr_insert_entries
   module procedure cuscr_insert_entries
   module procedure zuscr_insert_entries
end interface
 contains
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
      subroutine iuscr_insert_entries (a,val,indx,jndx,istat)
      implicit none
      integer  ,dimension(:),intent(in) ::val
      integer ,intent(in) ::a
      integer,intent(out)::istat
      integer,dimension(:),intent(in)::indx,jndx
      integer ::i
      istat=-1
      do i=1,size(val)
         call iuscr_insert_entry (a,val(i),indx(i),&
                                 jndx(i),istat)
         if(istat.ne.0) return
      end do
      end subroutine iuscr_insert_entries  
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
      subroutine suscr_insert_entries (a,val,indx,jndx,istat)
      implicit none
      real(KIND=sp)  ,dimension(:),intent(in) ::val
      integer ,intent(in) ::a
      integer,intent(out)::istat
      integer,dimension(:),intent(in)::indx,jndx
      integer ::i
      istat=-1
      do i=1,size(val)
         call suscr_insert_entry (a,val(i),indx(i),&
                                 jndx(i),istat)
         if(istat.ne.0) return
      end do
      end subroutine suscr_insert_entries  
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
      subroutine duscr_insert_entries (a,val,indx,jndx,istat)
      implicit none
      real(KIND=dp)  ,dimension(:),intent(in) ::val
      integer ,intent(in) ::a
      integer,intent(out)::istat
      integer,dimension(:),intent(in)::indx,jndx
      integer ::i
      istat=-1
      do i=1,size(val)
         call duscr_insert_entry (a,val(i),indx(i),&
                                 jndx(i),istat)
         if(istat.ne.0) return
      end do
      end subroutine duscr_insert_entries  
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
      subroutine cuscr_insert_entries (a,val,indx,jndx,istat)
      implicit none
      complex(KIND=sp)  ,dimension(:),intent(in) ::val
      integer ,intent(in) ::a
      integer,intent(out)::istat
      integer,dimension(:),intent(in)::indx,jndx
      integer ::i
      istat=-1
      do i=1,size(val)
         call cuscr_insert_entry (a,val(i),indx(i),&
                                 jndx(i),istat)
         if(istat.ne.0) return
      end do
      end subroutine cuscr_insert_entries  
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
      subroutine zuscr_insert_entries (a,val,indx,jndx,istat)
      implicit none
      complex(KIND=dp)  ,dimension(:),intent(in) ::val
      integer ,intent(in) ::a
      integer,intent(out)::istat
      integer,dimension(:),intent(in)::indx,jndx
      integer ::i
      istat=-1
      do i=1,size(val)
         call zuscr_insert_entry (a,val(i),indx(i),&
                                 jndx(i),istat)
         if(istat.ne.0) return
      end do
      end subroutine zuscr_insert_entries  
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
      end module mod_uscr_insert_entries
