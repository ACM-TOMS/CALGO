      module mod_usconv_csr2coo
      use properties
      use mod_conv_tools
      use representation_of_data
      contains
      subroutine usconv_csr2coo(a,ierr) 
      integer,intent(inout) :: a 
      integer ,intent(inout)::ierr
      integer :: res,s,rest
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
      if(isp_data %FIDA=='CSR') then
         call get_infoa(isp_data %INFOA ,'c',res,ierr)
         if(res.eq.COP_OF_SOURCE) then
            isp_data %FIDA='COO'
            s=size(isp_data %A)
            allocate(isp_data %IA2(s))
            isp_data %IA2= isp_data %IA1
            call PNTR_INV( isp_data %PE, isp_data %IA1)
            nullify(isp_data %PB)
            nullify(isp_data %PE)         
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
      if(ssp_data %FIDA=='CSR') then
         call get_infoa(ssp_data %INFOA ,'c',res,ierr)
         if(res.eq.COP_OF_SOURCE) then
            ssp_data %FIDA='COO'
            s=size(ssp_data %A)
            allocate(ssp_data %IA2(s))
            ssp_data %IA2= ssp_data %IA1
            call PNTR_INV( ssp_data %PE, ssp_data %IA1)
            nullify(ssp_data %PB)
            nullify(ssp_data %PE)         
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
      if(dsp_data %FIDA=='CSR') then
         call get_infoa(dsp_data %INFOA ,'c',res,ierr)
         if(res.eq.COP_OF_SOURCE) then
            dsp_data %FIDA='COO'
            s=size(dsp_data %A)
            allocate(dsp_data %IA2(s))
            dsp_data %IA2= dsp_data %IA1
            call PNTR_INV( dsp_data %PE, dsp_data %IA1)
            nullify(dsp_data %PB)
            nullify(dsp_data %PE)         
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
      if(csp_data %FIDA=='CSR') then
         call get_infoa(csp_data %INFOA ,'c',res,ierr)
         if(res.eq.COP_OF_SOURCE) then
            csp_data %FIDA='COO'
            s=size(csp_data %A)
            allocate(csp_data %IA2(s))
            csp_data %IA2= csp_data %IA1
            call PNTR_INV( csp_data %PE, csp_data %IA1)
            nullify(csp_data %PB)
            nullify(csp_data %PE)         
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
      if(zsp_data %FIDA=='CSR') then
         call get_infoa(zsp_data %INFOA ,'c',res,ierr)
         if(res.eq.COP_OF_SOURCE) then
            zsp_data %FIDA='COO'
            s=size(zsp_data %A)
            allocate(zsp_data %IA2(s))
            zsp_data %IA2= zsp_data %IA1
            call PNTR_INV( zsp_data %PE, zsp_data %IA1)
            nullify(zsp_data %PB)
            nullify(zsp_data %PE)         
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
      end subroutine usconv_csr2coo
      end module mod_usconv_csr2coo
