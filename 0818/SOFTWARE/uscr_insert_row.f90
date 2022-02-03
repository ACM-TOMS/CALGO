module mod_uscr_insert_row
use blas_sparse_namedconstants
use mod_uscr_insert_entry
interface uscr_insert_row
  module procedure iuscr_insert_row
  module procedure suscr_insert_row
  module procedure duscr_insert_row
  module procedure cuscr_insert_row
  module procedure zuscr_insert_row
end interface
 contains
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
      subroutine iuscr_insert_row  (a,i,val,jndx,istat)
      implicit none
      integer  ,dimension(:),intent(in) ::val
      integer ,intent(in) ::a,i
      integer,dimension(:),intent(in)::jndx
      integer ,intent(out)::istat
      integer ::k,s
      istat=-1
      s=size(val)
      do k=1,s
         call iuscr_insert_entry (a,val(k),i,jndx(k),istat)
         if (istat.ne.0) return
      end do
      end subroutine iuscr_insert_row 
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
      subroutine suscr_insert_row  (a,i,val,jndx,istat)
      implicit none
      real(KIND=sp)  ,dimension(:),intent(in) ::val
      integer ,intent(in) ::a,i
      integer,dimension(:),intent(in)::jndx
      integer ,intent(out)::istat
      integer ::k,s
      istat=-1
      s=size(val)
      do k=1,s
         call suscr_insert_entry (a,val(k),i,jndx(k),istat)
         if (istat.ne.0) return
      end do
      end subroutine suscr_insert_row 
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
      subroutine duscr_insert_row  (a,i,val,jndx,istat)
      implicit none
      real(KIND=dp)  ,dimension(:),intent(in) ::val
      integer ,intent(in) ::a,i
      integer,dimension(:),intent(in)::jndx
      integer ,intent(out)::istat
      integer ::k,s
      istat=-1
      s=size(val)
      do k=1,s
         call duscr_insert_entry (a,val(k),i,jndx(k),istat)
         if (istat.ne.0) return
      end do
      end subroutine duscr_insert_row 
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
      subroutine cuscr_insert_row  (a,i,val,jndx,istat)
      implicit none
      complex(KIND=sp)  ,dimension(:),intent(in) ::val
      integer ,intent(in) ::a,i
      integer,dimension(:),intent(in)::jndx
      integer ,intent(out)::istat
      integer ::k,s
      istat=-1
      s=size(val)
      do k=1,s
         call cuscr_insert_entry (a,val(k),i,jndx(k),istat)
         if (istat.ne.0) return
      end do
      end subroutine cuscr_insert_row 
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
      subroutine zuscr_insert_row  (a,i,val,jndx,istat)
      implicit none
      complex(KIND=dp)  ,dimension(:),intent(in) ::val
      integer ,intent(in) ::a,i
      integer,dimension(:),intent(in)::jndx
      integer ,intent(out)::istat
      integer ::k,s
      istat=-1
      s=size(val)
      do k=1,s
         call zuscr_insert_entry (a,val(k),i,jndx(k),istat)
         if (istat.ne.0) return
      end do
      end subroutine zuscr_insert_row 
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
end module mod_uscr_insert_row
