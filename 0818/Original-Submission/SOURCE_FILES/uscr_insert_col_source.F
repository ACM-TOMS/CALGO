      subroutine  XUSCR_INSERT_COL (a,j,val,indx,istat)
      implicit none
      DCOMPLEX ,dimension(:),intent(in) ::val
      integer ,intent(in) ::a,j
      integer ,intent(out) ::istat
      integer,dimension(:),intent(in)::indx
      integer ::i,s
      
      istat=-1
      s=size(val)
      do i=1,s
         call   XUSCR_INSERT_ENTRY(a,val(i),indx(i),j,istat)
         if(istat.ne.0) return
      end do
      
      end subroutine  XUSCR_INSERT_COL 
