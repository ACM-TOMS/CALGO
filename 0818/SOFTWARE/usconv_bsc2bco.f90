      module mod_usconv_bsc2bco
      use properties
      use mod_conv_tools
      use representation_of_data
      contains
      subroutine usconv_bsc2bco(a,ierr) 
      implicit none
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
      if(isp_data %FIDA=='BSC') then
         call get_infoa(isp_data %INFOA ,'c',res,ierr)
         if(res.eq.COP_OF_SOURCE) then
            isp_data %FIDA='BCO'
            s=size(isp_data %IA1)
            allocate(isp_data %IA2(s))
            call PNTR_INV( isp_data %PE, isp_data %IA2)
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
      if(ssp_data %FIDA=='BSC') then
         call get_infoa(ssp_data %INFOA ,'c',res,ierr)
         if(res.eq.COP_OF_SOURCE) then
            ssp_data %FIDA='BCO'
            s=size(ssp_data %IA1)
            allocate(ssp_data %IA2(s))
            call PNTR_INV( ssp_data %PE, ssp_data %IA2)
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
      if(dsp_data %FIDA=='BSC') then
         call get_infoa(dsp_data %INFOA ,'c',res,ierr)
         if(res.eq.COP_OF_SOURCE) then
            dsp_data %FIDA='BCO'
            s=size(dsp_data %IA1)
            allocate(dsp_data %IA2(s))
            call PNTR_INV( dsp_data %PE, dsp_data %IA2)
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
      if(csp_data %FIDA=='BSC') then
         call get_infoa(csp_data %INFOA ,'c',res,ierr)
         if(res.eq.COP_OF_SOURCE) then
            csp_data %FIDA='BCO'
            s=size(csp_data %IA1)
            allocate(csp_data %IA2(s))
            call PNTR_INV( csp_data %PE, csp_data %IA2)
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
      if(zsp_data %FIDA=='BSC') then
         call get_infoa(zsp_data %INFOA ,'c',res,ierr)
         if(res.eq.COP_OF_SOURCE) then
            zsp_data %FIDA='BCO'
            s=size(zsp_data %IA1)
            allocate(zsp_data %IA2(s))
            call PNTR_INV( zsp_data %PE, zsp_data %IA2)
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
      end subroutine usconv_bsc2bco
      end module mod_usconv_bsc2bco
