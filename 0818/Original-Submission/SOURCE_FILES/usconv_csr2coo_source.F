      call accessdata(XSP_DATA,a,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      if(XSP_DATA%FIDA=='CSR') then
         call get_infoa(XSP_DATA%INFOA ,'c',res,ierr)
         if(res.eq.COP_OF_SOURCE) then
            XSP_DATA%FIDA='COO'
            s=size(XSP_DATA%A)
            allocate(XSP_DATA%IA2(s))
            XSP_DATA%IA2=XSP_DATA%IA1
            call PNTR_INV( XSP_DATA%PE, XSP_DATA%IA1)
            nullify(XSP_DATA%PB)
            nullify(XSP_DATA%PE)         
         end if
      else 
         ierr = blas_error_param
         return
      end if
