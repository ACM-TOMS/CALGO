      module mod_uscr_end
      use mod_INSERTING
      use mod_INS_ROUTINER
      use properties
      contains
      subroutine uscr_end(a,istat)
      implicit none
      integer ,intent(inout)::a,istat
      integer::prpty,rest,b
      type(i_matrix),pointer ::ipmatrix
      type(d_matrix),pointer ::dpmatrix
      type(s_matrix),pointer ::spmatrix
      type(c_matrix),pointer ::cpmatrix
      type(z_matrix),pointer ::zpmatrix
      b=-a
      rest = modulo(b,no_of_types)
      select case(rest)
      case(ISP_MATRIX)
! **********************************************************************
      istat=-1
      call access_matrix(ipmatrix ,a,istat)  
      if(istat.ne.0) return
      prpty= ipmatrix %property  
      select case(ipmatrix %format)
      case('block')
         call iuscr_blockend (a,prpty,istat)
         if(istat.ne.0) return
      case('vblock')
         call iuscr_varend  (a,prpty,istat)
         if(istat.ne.0) return
      case('normal')
         call iuscr_normend (a,prpty,istat)       
         if(istat.ne.0) return
      case default
         istat = blas_error_param
         return
      end select 
!!*************************************************************************** 
! **********************************************************************
      case(SSP_MATRIX)
! **********************************************************************
      istat=-1
      call access_matrix(spmatrix ,a,istat)  
      if(istat.ne.0) return
      prpty= spmatrix %property  
      select case(spmatrix %format)
      case('block')
         call suscr_blockend (a,prpty,istat)
         if(istat.ne.0) return
      case('vblock')
         call suscr_varend  (a,prpty,istat)
         if(istat.ne.0) return
      case('normal')
         call suscr_normend (a,prpty,istat)       
         if(istat.ne.0) return
      case default
         istat = blas_error_param
         return
      end select 
! **********************************************************************
!!***************************************************************************
      case(DSP_MATRIX)
! **********************************************************************
      istat=-1
      call access_matrix(dpmatrix ,a,istat)  
      if(istat.ne.0) return
      prpty= dpmatrix %property  
      select case(dpmatrix %format)
      case('block')
         call duscr_blockend (a,prpty,istat)
         if(istat.ne.0) return
      case('vblock')
         call duscr_varend  (a,prpty,istat)
         if(istat.ne.0) return
      case('normal')
         call duscr_normend (a,prpty,istat)       
         if(istat.ne.0) return
      case default
         istat = blas_error_param
         return
      end select 
! **********************************************************************
!!***************************************************************************
      case(CSP_MATRIX)
! **********************************************************************
      istat=-1
      call access_matrix(cpmatrix ,a,istat)  
      if(istat.ne.0) return
      prpty= cpmatrix %property  
      select case(cpmatrix %format)
      case('block')
         call cuscr_blockend (a,prpty,istat)
         if(istat.ne.0) return
      case('vblock')
         call cuscr_varend  (a,prpty,istat)
         if(istat.ne.0) return
      case('normal')
         call cuscr_normend (a,prpty,istat)       
         if(istat.ne.0) return
      case default
         istat = blas_error_param
         return
      end select 
! **********************************************************************
!!*************************************************************************** 
      case(ZSP_MATRIX)
! **********************************************************************
      istat=-1
      call access_matrix(zpmatrix ,a,istat)  
      if(istat.ne.0) return
      prpty= zpmatrix %property  
      select case(zpmatrix %format)
      case('block')
         call zuscr_blockend (a,prpty,istat)
         if(istat.ne.0) return
      case('vblock')
         call zuscr_varend  (a,prpty,istat)
         if(istat.ne.0) return
      case('normal')
         call zuscr_normend (a,prpty,istat)       
         if(istat.ne.0) return
      case default
         istat = blas_error_param
         return
      end select 
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
      case default
         istat = blas_error_param
         return
      end select
      end subroutine uscr_end
      end module mod_uscr_end
