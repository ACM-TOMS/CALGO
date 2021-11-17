      module mod_usconv_bco2bsc
      use properties
      use mod_conv_tools
      use representation_of_data
      contains
      subroutine usconv_bco2bsc(a,ierr) 
      implicit none
      integer,intent(inout) :: a 
      integer ,intent(inout)::ierr
      integer :: res,col_dim_in_blocks,col_dim_of_block,rest
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
      if(isp_data %FIDA=='BCO') then
         call get_infoa(isp_data %INFOA ,'c',res,ierr)
         if(res.eq.COP_OF_SOURCE) then
            isp_data %FIDA='BSC'
            call get_infoa(isp_data %INFOA,'g',col_dim_in_blocks,ierr)
            call get_infoa(isp_data %INFOA,'e',col_dim_of_block,ierr)
            allocate(isp_data %PB(col_dim_in_blocks))
            allocate(isp_data %PE(col_dim_in_blocks))
            call ipre_usconv_bco2bsc (isp_data %A,isp_data %IA1,&
                       isp_data %IA2,col_dim_in_blocks,&
                       col_dim_of_block,isp_data %PB,isp_data %PE)
            nullify(isp_data %IA2)
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
      if(ssp_data %FIDA=='BCO') then
         call get_infoa(ssp_data %INFOA ,'c',res,ierr)
         if(res.eq.COP_OF_SOURCE) then
            ssp_data %FIDA='BSC'
            call get_infoa(ssp_data %INFOA,'g',col_dim_in_blocks,ierr)
            call get_infoa(ssp_data %INFOA,'e',col_dim_of_block,ierr)
            allocate(ssp_data %PB(col_dim_in_blocks))
            allocate(ssp_data %PE(col_dim_in_blocks))
            call spre_usconv_bco2bsc (ssp_data %A,ssp_data %IA1,&
                       ssp_data %IA2,col_dim_in_blocks,&
                       col_dim_of_block,ssp_data %PB,ssp_data %PE)
            nullify(ssp_data %IA2)
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
      if(dsp_data %FIDA=='BCO') then
         call get_infoa(dsp_data %INFOA ,'c',res,ierr)
         if(res.eq.COP_OF_SOURCE) then
            dsp_data %FIDA='BSC'
            call get_infoa(dsp_data %INFOA,'g',col_dim_in_blocks,ierr)
            call get_infoa(dsp_data %INFOA,'e',col_dim_of_block,ierr)
            allocate(dsp_data %PB(col_dim_in_blocks))
            allocate(dsp_data %PE(col_dim_in_blocks))
            call dpre_usconv_bco2bsc (dsp_data %A,dsp_data %IA1,&
                       dsp_data %IA2,col_dim_in_blocks,&
                       col_dim_of_block,dsp_data %PB,dsp_data %PE)
            nullify(dsp_data %IA2)
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
      if(csp_data %FIDA=='BCO') then
         call get_infoa(csp_data %INFOA ,'c',res,ierr)
         if(res.eq.COP_OF_SOURCE) then
            csp_data %FIDA='BSC'
            call get_infoa(csp_data %INFOA,'g',col_dim_in_blocks,ierr)
            call get_infoa(csp_data %INFOA,'e',col_dim_of_block,ierr)
            allocate(csp_data %PB(col_dim_in_blocks))
            allocate(csp_data %PE(col_dim_in_blocks))
            call cpre_usconv_bco2bsc (csp_data %A,csp_data %IA1,&
                       csp_data %IA2,col_dim_in_blocks,&
                       col_dim_of_block,csp_data %PB,csp_data %PE)
            nullify(csp_data %IA2)
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
      if(zsp_data %FIDA=='BCO') then
         call get_infoa(zsp_data %INFOA ,'c',res,ierr)
         if(res.eq.COP_OF_SOURCE) then
            zsp_data %FIDA='BSC'
            call get_infoa(zsp_data %INFOA,'g',col_dim_in_blocks,ierr)
            call get_infoa(zsp_data %INFOA,'e',col_dim_of_block,ierr)
            allocate(zsp_data %PB(col_dim_in_blocks))
            allocate(zsp_data %PE(col_dim_in_blocks))
            call zpre_usconv_bco2bsc (zsp_data %A,zsp_data %IA1,&
                       zsp_data %IA2,col_dim_in_blocks,&
                       col_dim_of_block,zsp_data %PB,zsp_data %PE)
            nullify(zsp_data %IA2)
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
      end subroutine usconv_bco2bsc
      end module mod_usconv_bco2bsc
