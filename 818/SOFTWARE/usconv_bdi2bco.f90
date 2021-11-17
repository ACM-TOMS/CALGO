      module mod_usconv_bdi2bco
      use properties
      use mod_conv_tools
      use representation_of_data
      contains
      subroutine usconv_bdi2bco(a,ierr) 
      implicit none
      integer,intent(inout) :: a 
      integer ,intent(inout)::ierr 
      integer :: res,BLDA,BNNZ,mb,kb,lb,rest
      type(ispmat),pointer :: isp_data
      type(sspmat),pointer :: ssp_data
      type(dspmat),pointer :: dsp_data
      type(cspmat),pointer :: csp_data
      type(zspmat),pointer :: zsp_data
      intrinsic floor
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
      if(isp_data %FIDA=='BDI') then
         call get_infoa(isp_data %INFOA ,'c',res,ierr)
         if(res.eq.COP_OF_SOURCE) then
            isp_data %FIDA='BCO'
            allocate(isp_data %IA2(2))
            call get_infoa(isp_data %INFOA,'d',lb,ierr) !row-dim of a block
            call get_infoa(isp_data %INFOA,'f',BLDA,ierr) !blocks per diagonal
            mb=floor(real(isp_data %M/lb))
            kb=floor(real(isp_data %K/lb))
            call set_infoa(isp_data %INFOA,'f',mb,ierr) !row-dim in blocks
            call set_infoa(isp_data %INFOA,'g',kb,ierr) !col-dim in blocks
            call ipre_usconv_bdi2bco (isp_data %A, isp_data %IA1,& 
                              isp_data %IA2,BLDA,BNNZ,lb)
            call set_infoa(isp_data %INFOA,'n',BNNZ,ierr)    
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
      if(ssp_data %FIDA=='BDI') then
         call get_infoa(ssp_data %INFOA ,'c',res,ierr)
         if(res.eq.COP_OF_SOURCE) then
            ssp_data %FIDA='BCO'
            allocate(ssp_data %IA2(2))
            call get_infoa(ssp_data %INFOA,'d',lb,ierr) !row-dim of a block
            call get_infoa(ssp_data %INFOA,'f',BLDA,ierr) !blocks per diagonal
            mb=floor(real(ssp_data %M/lb))
            kb=floor(real(ssp_data %K/lb))
            call set_infoa(ssp_data %INFOA,'f',mb,ierr) !row-dim in blocks
            call set_infoa(ssp_data %INFOA,'g',kb,ierr) !col-dim in blocks
            call spre_usconv_bdi2bco (ssp_data %A, ssp_data %IA1,& 
                              ssp_data %IA2,BLDA,BNNZ,lb)
            call set_infoa(ssp_data %INFOA,'n',BNNZ,ierr)    
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
      if(dsp_data %FIDA=='BDI') then
         call get_infoa(dsp_data %INFOA ,'c',res,ierr)
         if(res.eq.COP_OF_SOURCE) then
            dsp_data %FIDA='BCO'
            allocate(dsp_data %IA2(2))
            call get_infoa(dsp_data %INFOA,'d',lb,ierr) !row-dim of a block
            call get_infoa(dsp_data %INFOA,'f',BLDA,ierr) !blocks per diagonal
            mb=floor(real(dsp_data %M/lb))
            kb=floor(real(dsp_data %K/lb))
            call set_infoa(dsp_data %INFOA,'f',mb,ierr) !row-dim in blocks
            call set_infoa(dsp_data %INFOA,'g',kb,ierr) !col-dim in blocks
            call dpre_usconv_bdi2bco (dsp_data %A, dsp_data %IA1,& 
                              dsp_data %IA2,BLDA,BNNZ,lb)
            call set_infoa(dsp_data %INFOA,'n',BNNZ,ierr)    
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
      if(csp_data %FIDA=='BDI') then
         call get_infoa(csp_data %INFOA ,'c',res,ierr)
         if(res.eq.COP_OF_SOURCE) then
            csp_data %FIDA='BCO'
            allocate(csp_data %IA2(2))
            call get_infoa(csp_data %INFOA,'d',lb,ierr) !row-dim of a block
            call get_infoa(csp_data %INFOA,'f',BLDA,ierr) !blocks per diagonal
            mb=floor(real(csp_data %M/lb))
            kb=floor(real(csp_data %K/lb))
            call set_infoa(csp_data %INFOA,'f',mb,ierr) !row-dim in blocks
            call set_infoa(csp_data %INFOA,'g',kb,ierr) !col-dim in blocks
            call cpre_usconv_bdi2bco (csp_data %A, csp_data %IA1,& 
                              csp_data %IA2,BLDA,BNNZ,lb)
            call set_infoa(csp_data %INFOA,'n',BNNZ,ierr)    
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
      if(zsp_data %FIDA=='BDI') then
         call get_infoa(zsp_data %INFOA ,'c',res,ierr)
         if(res.eq.COP_OF_SOURCE) then
            zsp_data %FIDA='BCO'
            allocate(zsp_data %IA2(2))
            call get_infoa(zsp_data %INFOA,'d',lb,ierr) !row-dim of a block
            call get_infoa(zsp_data %INFOA,'f',BLDA,ierr) !blocks per diagonal
            mb=floor(real(zsp_data %M/lb))
            kb=floor(real(zsp_data %K/lb))
            call set_infoa(zsp_data %INFOA,'f',mb,ierr) !row-dim in blocks
            call set_infoa(zsp_data %INFOA,'g',kb,ierr) !col-dim in blocks
            call zpre_usconv_bdi2bco (zsp_data %A, zsp_data %IA1,& 
                              zsp_data %IA2,BLDA,BNNZ,lb)
            call set_infoa(zsp_data %INFOA,'n',BNNZ,ierr)    
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
      end subroutine usconv_bdi2bco
      end module mod_usconv_bdi2bco
