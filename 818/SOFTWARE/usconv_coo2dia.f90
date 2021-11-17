      module mod_usconv_coo2dia
      use properties
      use mod_conv_tools
      use representation_of_data
      contains
      subroutine usconv_coo2dia(a,ierr) 
      implicit none
      integer,intent(inout) :: a
      type(ispmat),pointer :: isp_data
      type(sspmat),pointer :: ssp_data
      type(dspmat),pointer :: dsp_data
      type(cspmat),pointer :: csp_data
      type(zspmat),pointer :: zsp_data
      integer ,intent(inout)::ierr
      integer :: res,LDA,NDIAG,nnz,rest
      rest = modulo(a,no_of_types)
      select case(rest)
      case(ISP_MATRIX)
! **********************************************************************
      ierr=-1
      call accessdata(isp_data ,a,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      if(isp_data %FIDA=='COO') then
         call get_infoa(isp_data %INFOA ,'c',res,ierr)
         if(res.eq.COP_OF_SOURCE) then
            isp_data %FIDA='DIA'
            call  ipre_usconv_coo2dia ( isp_data %M, isp_data %K, isp_data %A,& 
                 isp_data %IA1, isp_data %IA2,LDA,NDIAG)
            nullify(isp_data %IA2)
            nnz = count( isp_data %A.ne.0.)
            if(nnz.le.lda*ndiag*0.5) then
               ! Warning Many zeros stored !
            end if
            call set_infoa(isp_data %INFOA,'n',nnz,ierr)
            call set_infoa(isp_data %INFOA,'d',LDA,ierr) !row-dim of val
            call set_infoa(isp_data %INFOA,'e',NDIAG,ierr) !col-dim of val
         end if
      else 
         ierr = blas_error_param
         return
      end if
! **********************************************************************
      case(SSP_MATRIX)
! **********************************************************************
      ierr=-1
      call accessdata(ssp_data ,a,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      if(ssp_data %FIDA=='COO') then
         call get_infoa(ssp_data %INFOA ,'c',res,ierr)
         if(res.eq.COP_OF_SOURCE) then
            ssp_data %FIDA='DIA'
            call  spre_usconv_coo2dia ( ssp_data %M, ssp_data %K, ssp_data %A,& 
                 ssp_data %IA1, ssp_data %IA2,LDA,NDIAG)
            nullify(ssp_data %IA2)
            nnz = count( ssp_data %A.ne.0.)
            if(nnz.le.lda*ndiag*0.5) then
               ! Warning Many zeros stored !
            end if
            call set_infoa(ssp_data %INFOA,'n',nnz,ierr)
            call set_infoa(ssp_data %INFOA,'d',LDA,ierr) !row-dim of val
            call set_infoa(ssp_data %INFOA,'e',NDIAG,ierr) !col-dim of val
         end if
      else 
         ierr = blas_error_param
         return
      end if
! **********************************************************************
      case(DSP_MATRIX)
! **********************************************************************
      ierr=-1
      call accessdata(dsp_data ,a,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      if(dsp_data %FIDA=='COO') then
         call get_infoa(dsp_data %INFOA ,'c',res,ierr)
         if(res.eq.COP_OF_SOURCE) then
            dsp_data %FIDA='DIA'
            call  dpre_usconv_coo2dia ( dsp_data %M, dsp_data %K, dsp_data %A,& 
                 dsp_data %IA1, dsp_data %IA2,LDA,NDIAG)
            nullify(dsp_data %IA2)
            nnz = count( dsp_data %A.ne.0.)
            if(nnz.le.lda*ndiag*0.5) then
               ! Warning Many zeros stored !
            end if
            call set_infoa(dsp_data %INFOA,'n',nnz,ierr)
            call set_infoa(dsp_data %INFOA,'d',LDA,ierr) !row-dim of val
            call set_infoa(dsp_data %INFOA,'e',NDIAG,ierr) !col-dim of val
         end if
      else 
         ierr = blas_error_param
         return
      end if
! **********************************************************************
      case(CSP_MATRIX)
! **********************************************************************
      ierr=-1
      call accessdata(csp_data ,a,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      if(csp_data %FIDA=='COO') then
         call get_infoa(csp_data %INFOA ,'c',res,ierr)
         if(res.eq.COP_OF_SOURCE) then
            csp_data %FIDA='DIA'
            call  cpre_usconv_coo2dia ( csp_data %M, csp_data %K, csp_data %A,& 
                 csp_data %IA1, csp_data %IA2,LDA,NDIAG)
            nullify(csp_data %IA2)
            nnz = count( csp_data %A.ne.0.)
            if(nnz.le.lda*ndiag*0.5) then
               ! Warning Many zeros stored !
            end if
            call set_infoa(csp_data %INFOA,'n',nnz,ierr)
            call set_infoa(csp_data %INFOA,'d',LDA,ierr) !row-dim of val
            call set_infoa(csp_data %INFOA,'e',NDIAG,ierr) !col-dim of val
         end if
      else 
         ierr = blas_error_param
         return
      end if
! **********************************************************************
      case(ZSP_MATRIX)
! **********************************************************************
      ierr=-1
      call accessdata(zsp_data ,a,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      if(zsp_data %FIDA=='COO') then
         call get_infoa(zsp_data %INFOA ,'c',res,ierr)
         if(res.eq.COP_OF_SOURCE) then
            zsp_data %FIDA='DIA'
            call  zpre_usconv_coo2dia ( zsp_data %M, zsp_data %K, zsp_data %A,& 
                 zsp_data %IA1, zsp_data %IA2,LDA,NDIAG)
            nullify(zsp_data %IA2)
            nnz = count( zsp_data %A.ne.0.)
            if(nnz.le.lda*ndiag*0.5) then
               ! Warning Many zeros stored !
            end if
            call set_infoa(zsp_data %INFOA,'n',nnz,ierr)
            call set_infoa(zsp_data %INFOA,'d',LDA,ierr) !row-dim of val
            call set_infoa(zsp_data %INFOA,'e',NDIAG,ierr) !col-dim of val
         end if
      else 
         ierr = blas_error_param
         return
      end if
! **********************************************************************
      case default
         ierr = blas_error_param
         return
      end select
      end subroutine usconv_coo2dia
      end module mod_usconv_coo2dia
