         call accessdata(XSP_DATA,nmb,ierr)
         if (ierr.ne.0) then
            ierr = blas_error_param
            return
         end if
         call get_infoa(XSP_DATA%INFOA,'c',val,ierr)
         if (ierr.ne.0) then
            ierr = blas_error_param 
            return
         end if
         if(val.eq.COP_OF_SOURCE) then 
      ! *** Deallocate extra storage for copy of matrix *** !
         select case(XSP_DATA%FIDA)
         case('COO','BCO')
            deallocate(XSP_DATA%A,XSP_DATA%IA1,XSP_DATA%IA2,STAT=ierr)
         case('CSC','BSC')
            deallocate(XSP_DATA%A,XSP_DATA%IA1,XSP_DATA%pb,XSP_DATA%pe,&
                       STAT=ierr)
         case('CSR','BSR')
            deallocate(XSP_DATA%A,XSP_DATA%IA1,XSP_DATA%pb,XSP_DATA%pe,&
                       STAT=ierr)
         case('DIA','BDI')
            deallocate(XSP_DATA%A,XSP_DATA%IA1,STAT=ierr)
         case('VBR')
            deallocate(XSP_DATA%A,XSP_DATA%IA1,XSP_DATA%IA2,&
                   XSP_DATA%PB,XSP_DATA%PE,XSP_DATA%BP1,XSP_DATA%BP2&
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
      call DEL_XSP(nmb,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_memdeloc
         return
      end if
