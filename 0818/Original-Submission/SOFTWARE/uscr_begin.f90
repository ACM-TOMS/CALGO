module mod_uscr_begin
   use mod_INSERTING
   implicit none
 contains
! **********************************************************************
! **********************************************************************
      subroutine iuscr_begin (m,n,a,istat)
      implicit none
      integer ,intent(in) ::m,n
      integer ,intent(out)::a,istat
      integer ::nmb,mb
      type(i_matrix ),pointer :: ipmatrix 
      mb=1
      istat = -1
      if((m.le.0).or.(n.le.0)) then
         istat = blas_error_param
         return
      else
         call new_i_matrix (nmb,mb,istat)
         if (istat.ne.0) then
            istat = blas_error_memalloc
            return
         end if
         call  access_matrix(ipmatrix ,nmb,istat)
         if (istat.ne.0) then
            istat = blas_error_param
            return
      end if
      ipmatrix %DIM(1)=m !nb_of_rows
      ipmatrix %DIM(2)=n !nb_of_cols
      ipmatrix %format='normal'
      a=nmb
      end if
      end subroutine iuscr_begin 
! **********************************************************************
! **********************************************************************
      subroutine suscr_begin (m,n,a,istat)
      implicit none
      integer ,intent(in) ::m,n
      integer ,intent(out)::a,istat
      integer ::nmb,mb
      type(s_matrix ),pointer :: spmatrix 
      mb=1
      istat = -1
      if((m.le.0).or.(n.le.0)) then
         istat = blas_error_param
         return
      else
         call new_s_matrix (nmb,mb,istat)
         if (istat.ne.0) then
            istat = blas_error_memalloc
            return
         end if
         call  access_matrix(spmatrix ,nmb,istat)
         if (istat.ne.0) then
            istat = blas_error_param
            return
      end if
      spmatrix %DIM(1)=m !nb_of_rows
      spmatrix %DIM(2)=n !nb_of_cols
      spmatrix %format='normal'
      a=nmb
      end if
      end subroutine suscr_begin 
! **********************************************************************
! **********************************************************************
      subroutine duscr_begin (m,n,a,istat)
      implicit none
      integer ,intent(in) ::m,n
      integer ,intent(out)::a,istat
      integer ::nmb,mb
      type(d_matrix ),pointer :: dpmatrix 
      mb=1
      istat = -1
      if((m.le.0).or.(n.le.0)) then
         istat = blas_error_param
         return
      else
         call new_d_matrix (nmb,mb,istat)
         if (istat.ne.0) then
            istat = blas_error_memalloc
            return
         end if
         call  access_matrix(dpmatrix ,nmb,istat)
         if (istat.ne.0) then
            istat = blas_error_param
            return
      end if
      dpmatrix %DIM(1)=m !nb_of_rows
      dpmatrix %DIM(2)=n !nb_of_cols
      dpmatrix %format='normal'
      a=nmb
      end if
      end subroutine duscr_begin 
! **********************************************************************
! **********************************************************************
      subroutine cuscr_begin (m,n,a,istat)
      implicit none
      integer ,intent(in) ::m,n
      integer ,intent(out)::a,istat
      integer ::nmb,mb
      type(c_matrix ),pointer :: cpmatrix 
      mb=1
      istat = -1
      if((m.le.0).or.(n.le.0)) then
         istat = blas_error_param
         return
      else
         call new_c_matrix (nmb,mb,istat)
         if (istat.ne.0) then
            istat = blas_error_memalloc
            return
         end if
         call  access_matrix(cpmatrix ,nmb,istat)
         if (istat.ne.0) then
            istat = blas_error_param
            return
      end if
      cpmatrix %DIM(1)=m !nb_of_rows
      cpmatrix %DIM(2)=n !nb_of_cols
      cpmatrix %format='normal'
      a=nmb
      end if
      end subroutine cuscr_begin 
! **********************************************************************
! **********************************************************************
      subroutine zuscr_begin (m,n,a,istat)
      implicit none
      integer ,intent(in) ::m,n
      integer ,intent(out)::a,istat
      integer ::nmb,mb
      type(z_matrix ),pointer :: zpmatrix 
      mb=1
      istat = -1
      if((m.le.0).or.(n.le.0)) then
         istat = blas_error_param
         return
      else
         call new_z_matrix (nmb,mb,istat)
         if (istat.ne.0) then
            istat = blas_error_memalloc
            return
         end if
         call  access_matrix(zpmatrix ,nmb,istat)
         if (istat.ne.0) then
            istat = blas_error_param
            return
      end if
      zpmatrix %DIM(1)=m !nb_of_rows
      zpmatrix %DIM(2)=n !nb_of_cols
      zpmatrix %format='normal'
      a=nmb
      end if
      end subroutine zuscr_begin 
! **********************************************************************
! **********************************************************************
end module mod_uscr_begin
