      subroutine  XUSCR_INSERT_CLIQUE(a,val,indx,jndx,istat)
      implicit none
      DCOMPLEX ,dimension(:,:),intent(in) ::val
      integer ,intent(in) ::a
      integer ,intent(out) ::istat
      integer,dimension(:),intent(in)::indx,jndx
      integer ::i,j,s_row,s_col
      
      istat=-1
      s_row=size(indx)
      s_col=size(jndx)
      
      do j=1,s_col
         do i=1,s_row
            call   XUSCR_INSERT_ENTRY(a,val(i,j),&
                indx(i),jndx(j),istat)  
            if(istat.ne.0) return
         end do
      end do

      end subroutine  XUSCR_INSERT_CLIQUE
