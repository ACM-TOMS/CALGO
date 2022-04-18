      subroutine XUSCR_INSERT_ENTRIES(a,val,indx,jndx,istat)
      implicit none
      DCOMPLEX ,dimension(:),intent(in) ::val
      integer ,intent(in) ::a
      integer,intent(out)::istat
      integer,dimension(:),intent(in)::indx,jndx
      integer ::i

      istat=-1
      do i=1,size(val)
         call XUSCR_INSERT_ENTRY(a,val(i),indx(i),&
                                 jndx(i),istat)
         if(istat.ne.0) return
      end do

      end subroutine XUSCR_INSERT_ENTRIES 
