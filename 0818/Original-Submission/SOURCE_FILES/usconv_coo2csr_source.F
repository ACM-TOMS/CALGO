      call accessdata(XSP_DATA,a,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      if(XSP_DATA%FIDA=='COO') then
         call get_infoa(XSP_DATA%INFOA ,'c',res,ierr)
         if(res.eq.COP_OF_SOURCE) then
            XSP_DATA%FIDA='CSR'
            allocate(XSP_DATA%PB(XSP_DATA%K))
            allocate(XSP_DATA%PE(XSP_DATA%K))
            call XPRE_COO2CSR ( XSP_DATA%A, XSP_DATA%IA1, &
               XSP_DATA%IA2, XSP_DATA%M, XSP_DATA%PB, XSP_DATA%PE)
            nullify(XSP_DATA%IA2)                         
         end if
      else 
         ierr = blas_error_param
         return
      end if
