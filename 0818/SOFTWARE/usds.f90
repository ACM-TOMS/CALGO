      module mod_usds
! **********************************************************************
!     Author : C. Voemel
!     Date of last modification : 7.7.00
!     Description : RELEASES HANDLES, LOOKS FOR THE "FREEDOM OF MEMORY"
! **********************************************************************
      use representation_of_data
      use properties
      implicit none
      contains
      subroutine usds(nmb,ierr)
      implicit none
      intrinsic modulo
      integer, intent(in) :: nmb
      integer, intent(out) :: ierr
      type(ispmat),pointer :: isp_data
      type(sspmat),pointer :: ssp_data
      type(dspmat),pointer :: dsp_data
      type(cspmat),pointer :: csp_data
      type(zspmat),pointer :: zsp_data
      integer :: rest,val
      rest = modulo(nmb,no_of_types)
      select case(rest)
      case(ISP_MATRIX)
! **********************************************************************
         call accessdata(isp_data ,nmb,ierr)
         if (ierr.ne.0) then
            ierr = blas_error_param
            return
         end if
         call get_infoa(isp_data %INFOA,'c',val,ierr)
         if (ierr.ne.0) then
            ierr = blas_error_param 
            return
         end if
         if(val.eq.COP_OF_SOURCE) then 
      ! *** Deallocate extra storage for copy of matrix *** !
         select case(isp_data %FIDA)
         case('COO','BCO')
            deallocate(isp_data %A,isp_data %IA1,isp_data %IA2,STAT=ierr)
         case('CSC','BSC')
            deallocate(isp_data %A,isp_data %IA1,isp_data %pb,isp_data %pe,&
                       STAT=ierr)
         case('CSR','BSR')
            deallocate(isp_data %A,isp_data %IA1,isp_data %pb,isp_data %pe,&
                       STAT=ierr)
         case('DIA','BDI')
            deallocate(isp_data %A,isp_data %IA1,STAT=ierr)
         case('VBR')
            deallocate(isp_data %A,isp_data %IA1,isp_data %IA2,&
                   isp_data %PB,isp_data %PE,isp_data %BP1,isp_data %BP2&
                   ,STAT=ierr)
         case default
            ierr = blas_error_param
            return
         end select
         if(ierr.ne.0) then
            ierr = blas_error_memdeloc
            return
         end if
      end if
      call del_isp (nmb,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_memdeloc
         return
      end if
! **********************************************************************
      case(SSP_MATRIX)
! **********************************************************************
         call accessdata(ssp_data ,nmb,ierr)
         if (ierr.ne.0) then
            ierr = blas_error_param
            return
         end if
         call get_infoa(ssp_data %INFOA,'c',val,ierr)
         if (ierr.ne.0) then
            ierr = blas_error_param 
            return
         end if
         if(val.eq.COP_OF_SOURCE) then 
      ! *** Deallocate extra storage for copy of matrix *** !
         select case(ssp_data %FIDA)
         case('COO','BCO')
            deallocate(ssp_data %A,ssp_data %IA1,ssp_data %IA2,STAT=ierr)
         case('CSC','BSC')
            deallocate(ssp_data %A,ssp_data %IA1,ssp_data %pb,ssp_data %pe,&
                       STAT=ierr)
         case('CSR','BSR')
            deallocate(ssp_data %A,ssp_data %IA1,ssp_data %pb,ssp_data %pe,&
                       STAT=ierr)
         case('DIA','BDI')
            deallocate(ssp_data %A,ssp_data %IA1,STAT=ierr)
         case('VBR')
            deallocate(ssp_data %A,ssp_data %IA1,ssp_data %IA2,&
                   ssp_data %PB,ssp_data %PE,ssp_data %BP1,ssp_data %BP2&
                   ,STAT=ierr)
         case default
            ierr = blas_error_param
            return
         end select
         if(ierr.ne.0) then
            ierr = blas_error_memdeloc
            return
         end if
      end if
      call del_ssp (nmb,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_memdeloc
         return
      end if
! **********************************************************************
      case(DSP_MATRIX)
! **********************************************************************
         call accessdata(dsp_data ,nmb,ierr)
         if (ierr.ne.0) then
            ierr = blas_error_param
            return
         end if
         call get_infoa(dsp_data %INFOA,'c',val,ierr)
         if (ierr.ne.0) then
            ierr = blas_error_param 
            return
         end if
         if(val.eq.COP_OF_SOURCE) then 
      ! *** Deallocate extra storage for copy of matrix *** !
         select case(dsp_data %FIDA)
         case('COO','BCO')
            deallocate(dsp_data %A,dsp_data %IA1,dsp_data %IA2,STAT=ierr)
         case('CSC','BSC')
            deallocate(dsp_data %A,dsp_data %IA1,dsp_data %pb,dsp_data %pe,&
                       STAT=ierr)
         case('CSR','BSR')
            deallocate(dsp_data %A,dsp_data %IA1,dsp_data %pb,dsp_data %pe,&
                       STAT=ierr)
         case('DIA','BDI')
            deallocate(dsp_data %A,dsp_data %IA1,STAT=ierr)
         case('VBR')
            deallocate(dsp_data %A,dsp_data %IA1,dsp_data %IA2,&
                   dsp_data %PB,dsp_data %PE,dsp_data %BP1,dsp_data %BP2&
                   ,STAT=ierr)
         case default
            ierr = blas_error_param
            return
         end select
         if(ierr.ne.0) then
            ierr = blas_error_memdeloc
            return
         end if
      end if
      call del_dsp (nmb,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_memdeloc
         return
      end if
! **********************************************************************
      case(CSP_MATRIX)
! **********************************************************************
         call accessdata(csp_data ,nmb,ierr)
         if (ierr.ne.0) then
            ierr = blas_error_param
            return
         end if
         call get_infoa(csp_data %INFOA,'c',val,ierr)
         if (ierr.ne.0) then
            ierr = blas_error_param 
            return
         end if
         if(val.eq.COP_OF_SOURCE) then 
      ! *** Deallocate extra storage for copy of matrix *** !
         select case(csp_data %FIDA)
         case('COO','BCO')
            deallocate(csp_data %A,csp_data %IA1,csp_data %IA2,STAT=ierr)
         case('CSC','BSC')
            deallocate(csp_data %A,csp_data %IA1,csp_data %pb,csp_data %pe,&
                       STAT=ierr)
         case('CSR','BSR')
            deallocate(csp_data %A,csp_data %IA1,csp_data %pb,csp_data %pe,&
                       STAT=ierr)
         case('DIA','BDI')
            deallocate(csp_data %A,csp_data %IA1,STAT=ierr)
         case('VBR')
            deallocate(csp_data %A,csp_data %IA1,csp_data %IA2,&
                   csp_data %PB,csp_data %PE,csp_data %BP1,csp_data %BP2&
                   ,STAT=ierr)
         case default
            ierr = blas_error_param
            return
         end select
         if(ierr.ne.0) then
            ierr = blas_error_memdeloc
            return
         end if
      end if
      call del_csp (nmb,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_memdeloc
         return
      end if
! **********************************************************************
      case(ZSP_MATRIX)
! **********************************************************************
         call accessdata(zsp_data ,nmb,ierr)
         if (ierr.ne.0) then
            ierr = blas_error_param
            return
         end if
         call get_infoa(zsp_data %INFOA,'c',val,ierr)
         if (ierr.ne.0) then
            ierr = blas_error_param 
            return
         end if
         if(val.eq.COP_OF_SOURCE) then 
      ! *** Deallocate extra storage for copy of matrix *** !
         select case(zsp_data %FIDA)
         case('COO','BCO')
            deallocate(zsp_data %A,zsp_data %IA1,zsp_data %IA2,STAT=ierr)
         case('CSC','BSC')
            deallocate(zsp_data %A,zsp_data %IA1,zsp_data %pb,zsp_data %pe,&
                       STAT=ierr)
         case('CSR','BSR')
            deallocate(zsp_data %A,zsp_data %IA1,zsp_data %pb,zsp_data %pe,&
                       STAT=ierr)
         case('DIA','BDI')
            deallocate(zsp_data %A,zsp_data %IA1,STAT=ierr)
         case('VBR')
            deallocate(zsp_data %A,zsp_data %IA1,zsp_data %IA2,&
                   zsp_data %PB,zsp_data %PE,zsp_data %BP1,zsp_data %BP2&
                   ,STAT=ierr)
         case default
            ierr = blas_error_param
            return
         end select
         if(ierr.ne.0) then
            ierr = blas_error_memdeloc
            return
         end if
      end if
      call del_zsp (nmb,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_memdeloc
         return
      end if
! **********************************************************************
      case default
         ierr = blas_error_param 
      end select
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      end subroutine usds
      end module mod_usds
