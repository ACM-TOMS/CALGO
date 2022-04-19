      module mod_usconv_coo2csr
      use properties
      use mod_conv_tools
      use representation_of_data
      contains
      subroutine usconv_coo2csr(a,ierr) 
      implicit none
      integer,intent(inout) :: a 
      type(ispmat),pointer :: isp_data
      type(sspmat),pointer :: ssp_data
      type(dspmat),pointer :: dsp_data
      type(cspmat),pointer :: csp_data
      type(zspmat),pointer :: zsp_data
      integer ,intent(inout)::ierr
      integer :: res,rest
      ierr=-1
      rest = modulo(a,no_of_types)
      select case(rest)
      case(ISP_MATRIX) 
!!*************************************************************************** 
! **********************************************************************
      call accessdata(isp_data ,a,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      if(isp_data %FIDA=='COO') then
         call get_infoa(isp_data %INFOA ,'c',res,ierr)
         if(res.eq.COP_OF_SOURCE) then
            isp_data %FIDA='CSR'
            allocate(isp_data %PB(isp_data %K))
            allocate(isp_data %PE(isp_data %K))
            call ipre_usconv_coo2csr  ( isp_data %A, isp_data %IA1, &
               isp_data %IA2, isp_data %M, isp_data %PB, isp_data %PE)
            nullify(isp_data %IA2)                         
         end if
      else 
         ierr = blas_error_param
         return
      end if
      case(SSP_MATRIX) 
!!*************************************************************************** 
! **********************************************************************
      call accessdata(ssp_data ,a,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      if(ssp_data %FIDA=='COO') then
         call get_infoa(ssp_data %INFOA ,'c',res,ierr)
         if(res.eq.COP_OF_SOURCE) then
            ssp_data %FIDA='CSR'
            allocate(ssp_data %PB(ssp_data %K))
            allocate(ssp_data %PE(ssp_data %K))
            call spre_usconv_coo2csr  ( ssp_data %A, ssp_data %IA1, &
               ssp_data %IA2, ssp_data %M, ssp_data %PB, ssp_data %PE)
            nullify(ssp_data %IA2)                         
         end if
      else 
         ierr = blas_error_param
         return
      end if
! **********************************************************************
!!***************************************************************************
      case(DSP_MATRIX) 
! **********************************************************************
      call accessdata(dsp_data ,a,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      if(dsp_data %FIDA=='COO') then
         call get_infoa(dsp_data %INFOA ,'c',res,ierr)
         if(res.eq.COP_OF_SOURCE) then
            dsp_data %FIDA='CSR'
            allocate(dsp_data %PB(dsp_data %K))
            allocate(dsp_data %PE(dsp_data %K))
            call dpre_usconv_coo2csr  ( dsp_data %A, dsp_data %IA1, &
               dsp_data %IA2, dsp_data %M, dsp_data %PB, dsp_data %PE)
            nullify(dsp_data %IA2)                         
         end if
      else 
         ierr = blas_error_param
         return
      end if
! **********************************************************************
!!***************************************************************************
      case(CSP_MATRIX) 
! **********************************************************************
      call accessdata(csp_data ,a,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      if(csp_data %FIDA=='COO') then
         call get_infoa(csp_data %INFOA ,'c',res,ierr)
         if(res.eq.COP_OF_SOURCE) then
            csp_data %FIDA='CSR'
            allocate(csp_data %PB(csp_data %K))
            allocate(csp_data %PE(csp_data %K))
            call cpre_usconv_coo2csr  ( csp_data %A, csp_data %IA1, &
               csp_data %IA2, csp_data %M, csp_data %PB, csp_data %PE)
            nullify(csp_data %IA2)                         
         end if
      else 
         ierr = blas_error_param
         return
      end if
! **********************************************************************
!!***************************************************************************
      case(ZSP_MATRIX) 
! **********************************************************************
      call accessdata(zsp_data ,a,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      if(zsp_data %FIDA=='COO') then
         call get_infoa(zsp_data %INFOA ,'c',res,ierr)
         if(res.eq.COP_OF_SOURCE) then
            zsp_data %FIDA='CSR'
            allocate(zsp_data %PB(zsp_data %K))
            allocate(zsp_data %PE(zsp_data %K))
            call zpre_usconv_coo2csr  ( zsp_data %A, zsp_data %IA1, &
               zsp_data %IA2, zsp_data %M, zsp_data %PB, zsp_data %PE)
            nullify(zsp_data %IA2)                         
         end if
      else 
         ierr = blas_error_param
         return
      end if
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************  
      case default
         ierr = blas_error_param
         return
      end select
      end subroutine usconv_coo2csr
      end module mod_usconv_coo2csr
