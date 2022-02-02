      module mod_usconv_dia2coo
      use properties
      use mod_conv_tools
      use representation_of_data
      contains
      subroutine  usconv_dia2coo(a,ierr) 
      implicit none
      integer,intent(inout) :: a 
      integer ,intent(inout)::ierr
      integer :: res,LDA,NNZ,rest 
      type(ispmat),pointer :: isp_data
      type(sspmat),pointer :: ssp_data
      type(dspmat),pointer :: dsp_data
      type(cspmat),pointer :: csp_data
      type(zspmat),pointer :: zsp_data
      ierr=-1
      rest = modulo(a,no_of_types)
      select case(rest)
      case(ISP_MATRIX)
! **********************************************************************
      call accessdata(isp_data ,a,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      if(isp_data %FIDA=='DIA') then
         call get_infoa(isp_data %INFOA ,'c',res,ierr)
         if(res.eq.COP_OF_SOURCE) then
            isp_data %FIDA='COO'
            allocate(isp_data %IA2(2))
            call get_infoa(isp_data %INFOA ,'d',LDA,ierr) !row-dim of val
            call get_infoa(isp_data %INFOA ,'n',NNZ,ierr)
            call ipre_usconv_dia2coo (isp_data %A, isp_data %IA1, &
                 isp_data %IA2,LDA,NNZ)        
         end if
      else 
         ierr = blas_error_param
         return
      end if
! **********************************************************************
      case(SSP_MATRIX)
! **********************************************************************
      call accessdata(ssp_data ,a,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      if(ssp_data %FIDA=='DIA') then
         call get_infoa(ssp_data %INFOA ,'c',res,ierr)
         if(res.eq.COP_OF_SOURCE) then
            ssp_data %FIDA='COO'
            allocate(ssp_data %IA2(2))
            call get_infoa(ssp_data %INFOA ,'d',LDA,ierr) !row-dim of val
            call get_infoa(ssp_data %INFOA ,'n',NNZ,ierr)
            call spre_usconv_dia2coo (ssp_data %A, ssp_data %IA1, &
                 ssp_data %IA2,LDA,NNZ)        
         end if
      else 
         ierr = blas_error_param
         return
      end if
! **********************************************************************
      case(DSP_MATRIX)
! **********************************************************************
      call accessdata(dsp_data ,a,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      if(dsp_data %FIDA=='DIA') then
         call get_infoa(dsp_data %INFOA ,'c',res,ierr)
         if(res.eq.COP_OF_SOURCE) then
            dsp_data %FIDA='COO'
            allocate(dsp_data %IA2(2))
            call get_infoa(dsp_data %INFOA ,'d',LDA,ierr) !row-dim of val
            call get_infoa(dsp_data %INFOA ,'n',NNZ,ierr)
            call dpre_usconv_dia2coo (dsp_data %A, dsp_data %IA1, &
                 dsp_data %IA2,LDA,NNZ)        
         end if
      else 
         ierr = blas_error_param
         return
      end if
! **********************************************************************
      case(CSP_MATRIX)
! **********************************************************************
      call accessdata(csp_data ,a,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      if(csp_data %FIDA=='DIA') then
         call get_infoa(csp_data %INFOA ,'c',res,ierr)
         if(res.eq.COP_OF_SOURCE) then
            csp_data %FIDA='COO'
            allocate(csp_data %IA2(2))
            call get_infoa(csp_data %INFOA ,'d',LDA,ierr) !row-dim of val
            call get_infoa(csp_data %INFOA ,'n',NNZ,ierr)
            call cpre_usconv_dia2coo (csp_data %A, csp_data %IA1, &
                 csp_data %IA2,LDA,NNZ)        
         end if
      else 
         ierr = blas_error_param
         return
      end if
! **********************************************************************
      case(ZSP_MATRIX)
! **********************************************************************
      call accessdata(zsp_data ,a,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      if(zsp_data %FIDA=='DIA') then
         call get_infoa(zsp_data %INFOA ,'c',res,ierr)
         if(res.eq.COP_OF_SOURCE) then
            zsp_data %FIDA='COO'
            allocate(zsp_data %IA2(2))
            call get_infoa(zsp_data %INFOA ,'d',LDA,ierr) !row-dim of val
            call get_infoa(zsp_data %INFOA ,'n',NNZ,ierr)
            call zpre_usconv_dia2coo (zsp_data %A, zsp_data %IA1, &
                 zsp_data %IA2,LDA,NNZ)        
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
      end subroutine usconv_dia2coo
      end module mod_usconv_dia2coo
