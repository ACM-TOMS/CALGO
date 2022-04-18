module mod_uscr_variable_block_begin
  use mod_INSERTING
  use properties
  implicit none
 contains
! **********************************************************************
! **********************************************************************
      subroutine iuscr_variable_block_begin (Mb,Nb,k,l,a,istat)      
      implicit none
      integer ,intent(in) ::Mb,Nb
      integer,dimension(:),target,intent(in)::k,l
      integer ,intent(out)::a,istat
      integer ::nmb
      type(i_matrix ),pointer :: ipmatrix 
      istat = -1
      if((Mb.le.0).or.(Nb.le.0)) then
         istat = blas_error_param
         return
      else
         call  new_i_matrix (nmb,Mb, istat)
         if (istat.ne.0) return
         call  access_matrix(ipmatrix ,nmb, istat)
         if (istat.ne.0) return
         ipmatrix %DIM(3)=Mb     !nb_of_block_rows
         ipmatrix %DIM(4)=Nb     !nb_of_block_cols
         ipmatrix %sub_rows=>k
         ipmatrix %sub_cols=>l
         ipmatrix %trb=1
         ipmatrix %tre=1
         ipmatrix %format='vblock'
         a=nmb
      end if
      istat = 0
      end subroutine iuscr_variable_block_begin 
! **********************************************************************
! **********************************************************************
      subroutine suscr_variable_block_begin (Mb,Nb,k,l,a,istat)      
      implicit none
      integer ,intent(in) ::Mb,Nb
      integer,dimension(:),target,intent(in)::k,l
      integer ,intent(out)::a,istat
      integer ::nmb
      type(s_matrix ),pointer :: spmatrix 
      istat = -1
      if((Mb.le.0).or.(Nb.le.0)) then
         istat = blas_error_param
         return
      else
         call  new_s_matrix (nmb,Mb, istat)
         if (istat.ne.0) return
         call  access_matrix(spmatrix ,nmb, istat)
         if (istat.ne.0) return
         spmatrix %DIM(3)=Mb     !nb_of_block_rows
         spmatrix %DIM(4)=Nb     !nb_of_block_cols
         spmatrix %sub_rows=>k
         spmatrix %sub_cols=>l
         spmatrix %trb=1
         spmatrix %tre=1
         spmatrix %format='vblock'
         a=nmb
      end if
      istat = 0
      end subroutine suscr_variable_block_begin 
! **********************************************************************
! **********************************************************************
      subroutine duscr_variable_block_begin (Mb,Nb,k,l,a,istat)      
      implicit none
      integer ,intent(in) ::Mb,Nb
      integer,dimension(:),target,intent(in)::k,l
      integer ,intent(out)::a,istat
      integer ::nmb
      type(d_matrix ),pointer :: dpmatrix 
      istat = -1
      if((Mb.le.0).or.(Nb.le.0)) then
         istat = blas_error_param
         return
      else
         call  new_d_matrix (nmb,Mb, istat)
         if (istat.ne.0) return
         call  access_matrix(dpmatrix ,nmb, istat)
         if (istat.ne.0) return
         dpmatrix %DIM(3)=Mb     !nb_of_block_rows
         dpmatrix %DIM(4)=Nb     !nb_of_block_cols
         dpmatrix %sub_rows=>k
         dpmatrix %sub_cols=>l
         dpmatrix %trb=1
         dpmatrix %tre=1
         dpmatrix %format='vblock'
         a=nmb
      end if
      istat = 0
      end subroutine duscr_variable_block_begin 
! **********************************************************************
! **********************************************************************
      subroutine cuscr_variable_block_begin (Mb,Nb,k,l,a,istat)      
      implicit none
      integer ,intent(in) ::Mb,Nb
      integer,dimension(:),target,intent(in)::k,l
      integer ,intent(out)::a,istat
      integer ::nmb
      type(c_matrix ),pointer :: cpmatrix 
      istat = -1
      if((Mb.le.0).or.(Nb.le.0)) then
         istat = blas_error_param
         return
      else
         call  new_c_matrix (nmb,Mb, istat)
         if (istat.ne.0) return
         call  access_matrix(cpmatrix ,nmb, istat)
         if (istat.ne.0) return
         cpmatrix %DIM(3)=Mb     !nb_of_block_rows
         cpmatrix %DIM(4)=Nb     !nb_of_block_cols
         cpmatrix %sub_rows=>k
         cpmatrix %sub_cols=>l
         cpmatrix %trb=1
         cpmatrix %tre=1
         cpmatrix %format='vblock'
         a=nmb
      end if
      istat = 0
      end subroutine cuscr_variable_block_begin 
! **********************************************************************
! **********************************************************************
      subroutine zuscr_variable_block_begin (Mb,Nb,k,l,a,istat)      
      implicit none
      integer ,intent(in) ::Mb,Nb
      integer,dimension(:),target,intent(in)::k,l
      integer ,intent(out)::a,istat
      integer ::nmb
      type(z_matrix ),pointer :: zpmatrix 
      istat = -1
      if((Mb.le.0).or.(Nb.le.0)) then
         istat = blas_error_param
         return
      else
         call  new_z_matrix (nmb,Mb, istat)
         if (istat.ne.0) return
         call  access_matrix(zpmatrix ,nmb, istat)
         if (istat.ne.0) return
         zpmatrix %DIM(3)=Mb     !nb_of_block_rows
         zpmatrix %DIM(4)=Nb     !nb_of_block_cols
         zpmatrix %sub_rows=>k
         zpmatrix %sub_cols=>l
         zpmatrix %trb=1
         zpmatrix %tre=1
         zpmatrix %format='vblock'
         a=nmb
      end if
      istat = 0
      end subroutine zuscr_variable_block_begin 
! **********************************************************************
! **********************************************************************
end module mod_uscr_variable_block_begin
