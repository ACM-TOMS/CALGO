      subroutine XUSCR_BEGIN(m,n,a,istat)
      implicit none
      integer ,intent(in) ::m,n
      integer ,intent(out)::a,istat
      integer ::nmb,mb
      type(X_MATRIX),pointer :: XPMATRIX
      
      mb=1
      istat = -1
      if((m.le.0).or.(n.le.0)) then
         istat = blas_error_param
         return
      else
         call NEW_X_MATRIX(nmb,mb,istat)
         if (istat.ne.0) then
            istat = blas_error_memalloc
            return
         end if
      
         call  access_matrix(XPMATRIX,nmb,istat)
         if (istat.ne.0) then
            istat = blas_error_param
            return
      end if
  
      XPMATRIX%DIM(1)=m !nb_of_rows
      XPMATRIX%DIM(2)=n !nb_of_cols
      XPMATRIX%format='normal'
      a=nmb
      end if

      end subroutine XUSCR_BEGIN
