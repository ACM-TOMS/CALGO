      module mod_usconv_bco2bdi
      use properties
      use mod_conv_tools
      use representation_of_data
      contains
      subroutine usconv_bco2bdi(a,ierr) 
      implicit none
      integer,intent(inout) :: a 
      integer ,intent(inout)::ierr
      integer :: res,BLDA,BNDIAG,mb,kb,lb,rest   
      type(ispmat),pointer :: isp_data
      type(sspmat),pointer :: ssp_data
      type(dspmat),pointer :: dsp_data
      type(cspmat),pointer :: csp_data
      type(zspmat),pointer :: zsp_data
      intrinsic min   
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
            isp_data %FIDA='BDI'
            call get_infoa(isp_data %INFOA ,'e',lb,ierr)
            call get_infoa(isp_data %INFOA ,'f',mb,ierr)
            call get_infoa(isp_data %INFOA ,'g',kb,ierr)
            BLDA=min(mb,kb)
            call ipre_usconv_bco2bdi (mb,kb,lb,isp_data %A,&
                       isp_data %IA1,isp_data %IA2,BLDA,BNDIAG)     
            nullify(isp_data %IA2)
            call set_infoa(isp_data %INFOA,'d',lb,ierr) !row-dim of a block
            call set_infoa(isp_data %INFOA,'e',lb,ierr) !col-dim of a block
            call set_infoa(isp_data %INFOA,'f',BLDA,ierr) !blocks per diagonal
            call set_infoa(isp_data %INFOA,'g',BNDIAG,ierr) !no of diagonals
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
            ssp_data %FIDA='BDI'
            call get_infoa(ssp_data %INFOA ,'e',lb,ierr)
            call get_infoa(ssp_data %INFOA ,'f',mb,ierr)
            call get_infoa(ssp_data %INFOA ,'g',kb,ierr)
            BLDA=min(mb,kb)
            call spre_usconv_bco2bdi (mb,kb,lb,ssp_data %A,&
                       ssp_data %IA1,ssp_data %IA2,BLDA,BNDIAG)     
            nullify(ssp_data %IA2)
            call set_infoa(ssp_data %INFOA,'d',lb,ierr) !row-dim of a block
            call set_infoa(ssp_data %INFOA,'e',lb,ierr) !col-dim of a block
            call set_infoa(ssp_data %INFOA,'f',BLDA,ierr) !blocks per diagonal
            call set_infoa(ssp_data %INFOA,'g',BNDIAG,ierr) !no of diagonals
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
            dsp_data %FIDA='BDI'
            call get_infoa(dsp_data %INFOA ,'e',lb,ierr)
            call get_infoa(dsp_data %INFOA ,'f',mb,ierr)
            call get_infoa(dsp_data %INFOA ,'g',kb,ierr)
            BLDA=min(mb,kb)
            call dpre_usconv_bco2bdi (mb,kb,lb,dsp_data %A,&
                       dsp_data %IA1,dsp_data %IA2,BLDA,BNDIAG)     
            nullify(dsp_data %IA2)
            call set_infoa(dsp_data %INFOA,'d',lb,ierr) !row-dim of a block
            call set_infoa(dsp_data %INFOA,'e',lb,ierr) !col-dim of a block
            call set_infoa(dsp_data %INFOA,'f',BLDA,ierr) !blocks per diagonal
            call set_infoa(dsp_data %INFOA,'g',BNDIAG,ierr) !no of diagonals
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
            csp_data %FIDA='BDI'
            call get_infoa(csp_data %INFOA ,'e',lb,ierr)
            call get_infoa(csp_data %INFOA ,'f',mb,ierr)
            call get_infoa(csp_data %INFOA ,'g',kb,ierr)
            BLDA=min(mb,kb)
            call cpre_usconv_bco2bdi (mb,kb,lb,csp_data %A,&
                       csp_data %IA1,csp_data %IA2,BLDA,BNDIAG)     
            nullify(csp_data %IA2)
            call set_infoa(csp_data %INFOA,'d',lb,ierr) !row-dim of a block
            call set_infoa(csp_data %INFOA,'e',lb,ierr) !col-dim of a block
            call set_infoa(csp_data %INFOA,'f',BLDA,ierr) !blocks per diagonal
            call set_infoa(csp_data %INFOA,'g',BNDIAG,ierr) !no of diagonals
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
            zsp_data %FIDA='BDI'
            call get_infoa(zsp_data %INFOA ,'e',lb,ierr)
            call get_infoa(zsp_data %INFOA ,'f',mb,ierr)
            call get_infoa(zsp_data %INFOA ,'g',kb,ierr)
            BLDA=min(mb,kb)
            call zpre_usconv_bco2bdi (mb,kb,lb,zsp_data %A,&
                       zsp_data %IA1,zsp_data %IA2,BLDA,BNDIAG)     
            nullify(zsp_data %IA2)
            call set_infoa(zsp_data %INFOA,'d',lb,ierr) !row-dim of a block
            call set_infoa(zsp_data %INFOA,'e',lb,ierr) !col-dim of a block
            call set_infoa(zsp_data %INFOA,'f',BLDA,ierr) !blocks per diagonal
            call set_infoa(zsp_data %INFOA,'g',BNDIAG,ierr) !no of diagonals
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
      end subroutine usconv_bco2bdi
      end module mod_usconv_bco2bdi
