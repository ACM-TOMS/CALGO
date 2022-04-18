      subroutine XUSCR_INSERT_ENTRY (a,val,i,j,istat)
      implicit none
      DCOMPLEX ,intent(in) ::val
      integer ,intent(in) ::a,i,j
      integer,intent(out)::istat
      type(X_MATRIX),pointer ::pmatrix
  
      istat=-1
      call access_matrix(pmatrix,a,istat) 
      if (istat.ne.0) return
   
      select case(pmatrix%format)
      case('block')
         call XINS_bl_entr(pmatrix,val,i,j,istat)
      case('vblock')
         call  XINS_varbl_entr(pmatrix,val,i,j,istat) 
      case('normal')
         call XINS_entry(pmatrix,val,i,j,istat)
      case default
         istat = blas_error_param
         return
      end select

      end subroutine XUSCR_INSERT_ENTRY 
