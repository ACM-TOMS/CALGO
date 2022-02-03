      subroutine XUSCR_V_B_BEGIN(Mb,Nb,k,l,a,istat)      
      implicit none
      integer ,intent(in) ::Mb,Nb
      integer,dimension(:),target,intent(in)::k,l
      integer ,intent(out)::a,istat
      integer ::nmb
      type(X_MATRIX),pointer :: XPMATRIX
      
      istat = -1
      if((Mb.le.0).or.(Nb.le.0)) then
         istat = blas_error_param
         return
      else
         call  NEW_X_MATRIX(nmb,Mb, istat)
         if (istat.ne.0) return
         call  access_matrix(XPMATRIX,nmb, istat)
         if (istat.ne.0) return
         XPMATRIX%DIM(3)=Mb     !nb_of_block_rows
         XPMATRIX%DIM(4)=Nb     !nb_of_block_cols
         XPMATRIX%sub_rows=>k
         XPMATRIX%sub_cols=>l
         XPMATRIX%trb=1
         XPMATRIX%tre=1
         XPMATRIX%format='vblock'
         a=nmb
      end if
      istat = 0
      end subroutine XUSCR_V_B_BEGIN
