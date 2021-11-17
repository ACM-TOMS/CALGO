      subroutine  XUSCR_BLOCK_BEGIN(Mb,Nb,k,l,a,istat)
      implicit none
      integer ,intent(in) ::Mb,Nb,k,l
      integer ,intent(out)::a,istat
      integer ::nmb,m
      type(X_MATRIX),pointer :: XPMATRIX
      
      m=1
      istat = -1
      if((Mb.le.0).or.(Nb.le.0)) then
         istat = blas_error_param
         return
      else      
         call  NEW_X_MATRIX(nmb,m,istat)
         if (istat.ne.0) return      
         call  access_matrix(XPMATRIX,nmb, istat)
         if (istat.ne.0) return
         XPMATRIX%DIM(3)=Mb     !nb_of_block_rows
         XPMATRIX%DIM(4)=Nb     !nb_of_block_cols
         XPMATRIX%DIM(5)=k      !nb_of_rows_in_block
         XPMATRIX%DIM(6)=l      !nb_of_cols_in_block
         XPMATRIX%format='block'
         a=nmb
      end if
      istat = 0
      end subroutine XUSCR_BLOCK_BEGIN
