      subroutine XUSCR_INSERT_ROW (a,i,val,jndx,istat)
      implicit none
      DCOMPLEX ,dimension(:),intent(in) ::val
      integer ,intent(in) ::a,i
      integer,dimension(:),intent(in)::jndx
      integer ,intent(out)::istat
      integer ::k,s
      
      istat=-1
      s=size(val)
      do k=1,s
         call XUSCR_INSERT_ENTRY(a,val(k),i,jndx(k),istat)
         if (istat.ne.0) return
      end do
      
      end subroutine XUSCR_INSERT_ROW
