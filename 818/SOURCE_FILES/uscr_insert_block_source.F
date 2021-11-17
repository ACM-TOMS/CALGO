      subroutine  XUSCR_INSERT_BLOCK(a,val,bi,bj,istat)
      implicit none
      DCOMPLEX ,dimension(:,:),intent(in) ::val
      integer ,intent(in) ::a,bi,bj
      integer,intent(out)::istat
      type(X_MATRIX),pointer ::pmatrix

      istat=-1
      call access_matrix(pmatrix,a,istat) 
      if(istat.ne.0) return

      select case(pmatrix%format)
      case('block')
         call XINS_block(pmatrix,val,bi,bj,istat)
      case('vblock')
         call  XINS_varblock(pmatrix,val,bi,bj,istat)
      case default
         istat = blas_error_param    
         return
      end select
  
      end subroutine  XUSCR_INSERT_BLOCK
