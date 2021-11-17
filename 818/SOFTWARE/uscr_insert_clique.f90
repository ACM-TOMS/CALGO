module mod_uscr_insert_clique
use blas_sparse_namedconstants
use mod_uscr_insert_entry
interface uscr_insert_clique
module procedure iuscr_insert_clique
module procedure suscr_insert_clique
module procedure duscr_insert_clique
module procedure cuscr_insert_clique
module procedure zuscr_insert_clique
end interface
contains
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
      subroutine  iuscr_insert_clique (a,val,indx,jndx,istat)
      implicit none
      integer  ,dimension(:,:),intent(in) ::val
      integer ,intent(in) ::a
      integer ,intent(out) ::istat
      integer,dimension(:),intent(in)::indx,jndx
      integer ::i,j,s_row,s_col
      istat=-1
      s_row=size(indx)
      s_col=size(jndx)
      do j=1,s_col
         do i=1,s_row
            call   iuscr_insert_entry (a,val(i,j),&
                indx(i),jndx(j),istat)  
            if(istat.ne.0) return
         end do
      end do
      end subroutine  iuscr_insert_clique 
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
      subroutine  suscr_insert_clique (a,val,indx,jndx,istat)
      implicit none
      real(KIND=sp)  ,dimension(:,:),intent(in) ::val
      integer ,intent(in) ::a
      integer ,intent(out) ::istat
      integer,dimension(:),intent(in)::indx,jndx
      integer ::i,j,s_row,s_col
      istat=-1
      s_row=size(indx)
      s_col=size(jndx)
      do j=1,s_col
         do i=1,s_row
            call   suscr_insert_entry (a,val(i,j),&
                indx(i),jndx(j),istat)  
            if(istat.ne.0) return
         end do
      end do
      end subroutine  suscr_insert_clique 
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
      subroutine  duscr_insert_clique (a,val,indx,jndx,istat)
      implicit none
      real(KIND=dp)  ,dimension(:,:),intent(in) ::val
      integer ,intent(in) ::a
      integer ,intent(out) ::istat
      integer,dimension(:),intent(in)::indx,jndx
      integer ::i,j,s_row,s_col
      istat=-1
      s_row=size(indx)
      s_col=size(jndx)
      do j=1,s_col
         do i=1,s_row
            call   duscr_insert_entry (a,val(i,j),&
                indx(i),jndx(j),istat)  
            if(istat.ne.0) return
         end do
      end do
      end subroutine  duscr_insert_clique 
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
      subroutine  cuscr_insert_clique (a,val,indx,jndx,istat)
      implicit none
      complex(KIND=sp)  ,dimension(:,:),intent(in) ::val
      integer ,intent(in) ::a
      integer ,intent(out) ::istat
      integer,dimension(:),intent(in)::indx,jndx
      integer ::i,j,s_row,s_col
      istat=-1
      s_row=size(indx)
      s_col=size(jndx)
      do j=1,s_col
         do i=1,s_row
            call   cuscr_insert_entry (a,val(i,j),&
                indx(i),jndx(j),istat)  
            if(istat.ne.0) return
         end do
      end do
      end subroutine  cuscr_insert_clique 
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
      subroutine  zuscr_insert_clique (a,val,indx,jndx,istat)
      implicit none
      complex(KIND=dp)  ,dimension(:,:),intent(in) ::val
      integer ,intent(in) ::a
      integer ,intent(out) ::istat
      integer,dimension(:),intent(in)::indx,jndx
      integer ::i,j,s_row,s_col
      istat=-1
      s_row=size(indx)
      s_col=size(jndx)
      do j=1,s_col
         do i=1,s_row
            call   zuscr_insert_entry (a,val(i,j),&
                indx(i),jndx(j),istat)  
            if(istat.ne.0) return
         end do
      end do
      end subroutine  zuscr_insert_clique 
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
end module mod_uscr_insert_clique
