module mod_uscr_block_begin
  use mod_INSERTING
  use properties    
  implicit none
 contains
! **********************************************************************
! **********************************************************************
      subroutine  iuscr_block_begin (Mb,Nb,k,l,a,istat)
      implicit none
      integer ,intent(in) ::Mb,Nb,k,l
      integer ,intent(out)::a,istat
      integer ::nmb,m
      type(i_matrix ),pointer :: ipmatrix 
      m=1
      istat = -1
      if((Mb.le.0).or.(Nb.le.0)) then
         istat = blas_error_param
         return
      else      
         call  new_i_matrix (nmb,m,istat)
         if (istat.ne.0) return      
         call  access_matrix(ipmatrix ,nmb, istat)
         if (istat.ne.0) return
         ipmatrix %DIM(3)=Mb     !nb_of_block_rows
         ipmatrix %DIM(4)=Nb     !nb_of_block_cols
         ipmatrix %DIM(5)=k      !nb_of_rows_in_block
         ipmatrix %DIM(6)=l      !nb_of_cols_in_block
         ipmatrix %format='block'
         a=nmb
      end if
      istat = 0
      end subroutine iuscr_block_begin 
! **********************************************************************
! **********************************************************************
      subroutine  suscr_block_begin (Mb,Nb,k,l,a,istat)
      implicit none
      integer ,intent(in) ::Mb,Nb,k,l
      integer ,intent(out)::a,istat
      integer ::nmb,m
      type(s_matrix ),pointer :: spmatrix 
      m=1
      istat = -1
      if((Mb.le.0).or.(Nb.le.0)) then
         istat = blas_error_param
         return
      else      
         call  new_s_matrix (nmb,m,istat)
         if (istat.ne.0) return      
         call  access_matrix(spmatrix ,nmb, istat)
         if (istat.ne.0) return
         spmatrix %DIM(3)=Mb     !nb_of_block_rows
         spmatrix %DIM(4)=Nb     !nb_of_block_cols
         spmatrix %DIM(5)=k      !nb_of_rows_in_block
         spmatrix %DIM(6)=l      !nb_of_cols_in_block
         spmatrix %format='block'
         a=nmb
      end if
      istat = 0
      end subroutine suscr_block_begin 
! **********************************************************************
! **********************************************************************
      subroutine  duscr_block_begin (Mb,Nb,k,l,a,istat)
      implicit none
      integer ,intent(in) ::Mb,Nb,k,l
      integer ,intent(out)::a,istat
      integer ::nmb,m
      type(d_matrix ),pointer :: dpmatrix 
      m=1
      istat = -1
      if((Mb.le.0).or.(Nb.le.0)) then
         istat = blas_error_param
         return
      else      
         call  new_d_matrix (nmb,m,istat)
         if (istat.ne.0) return      
         call  access_matrix(dpmatrix ,nmb, istat)
         if (istat.ne.0) return
         dpmatrix %DIM(3)=Mb     !nb_of_block_rows
         dpmatrix %DIM(4)=Nb     !nb_of_block_cols
         dpmatrix %DIM(5)=k      !nb_of_rows_in_block
         dpmatrix %DIM(6)=l      !nb_of_cols_in_block
         dpmatrix %format='block'
         a=nmb
      end if
      istat = 0
      end subroutine duscr_block_begin 
! **********************************************************************
! **********************************************************************
      subroutine  cuscr_block_begin (Mb,Nb,k,l,a,istat)
      implicit none
      integer ,intent(in) ::Mb,Nb,k,l
      integer ,intent(out)::a,istat
      integer ::nmb,m
      type(c_matrix ),pointer :: cpmatrix 
      m=1
      istat = -1
      if((Mb.le.0).or.(Nb.le.0)) then
         istat = blas_error_param
         return
      else      
         call  new_c_matrix (nmb,m,istat)
         if (istat.ne.0) return      
         call  access_matrix(cpmatrix ,nmb, istat)
         if (istat.ne.0) return
         cpmatrix %DIM(3)=Mb     !nb_of_block_rows
         cpmatrix %DIM(4)=Nb     !nb_of_block_cols
         cpmatrix %DIM(5)=k      !nb_of_rows_in_block
         cpmatrix %DIM(6)=l      !nb_of_cols_in_block
         cpmatrix %format='block'
         a=nmb
      end if
      istat = 0
      end subroutine cuscr_block_begin 
! **********************************************************************
! **********************************************************************
      subroutine  zuscr_block_begin (Mb,Nb,k,l,a,istat)
      implicit none
      integer ,intent(in) ::Mb,Nb,k,l
      integer ,intent(out)::a,istat
      integer ::nmb,m
      type(z_matrix ),pointer :: zpmatrix 
      m=1
      istat = -1
      if((Mb.le.0).or.(Nb.le.0)) then
         istat = blas_error_param
         return
      else      
         call  new_z_matrix (nmb,m,istat)
         if (istat.ne.0) return      
         call  access_matrix(zpmatrix ,nmb, istat)
         if (istat.ne.0) return
         zpmatrix %DIM(3)=Mb     !nb_of_block_rows
         zpmatrix %DIM(4)=Nb     !nb_of_block_cols
         zpmatrix %DIM(5)=k      !nb_of_rows_in_block
         zpmatrix %DIM(6)=l      !nb_of_cols_in_block
         zpmatrix %format='block'
         a=nmb
      end if
      istat = 0
      end subroutine zuscr_block_begin 
! **********************************************************************
! **********************************************************************
end module mod_uscr_block_begin
