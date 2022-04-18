      call accessdata(XSP_DATA,a,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      if(XSP_DATA%FIDA=='BSC') then
         call get_infoa(XSP_DATA%INFOA ,'c',res,ierr)
         if(res.eq.COP_OF_SOURCE) then
            XSP_DATA%FIDA='BCO'
            s=size(XSP_DATA%IA1)
            allocate(XSP_DATA%IA2(s))
            call PNTR_INV( XSP_DATA%PE, XSP_DATA%IA2)
            nullify(XSP_DATA%PB)
            nullify(XSP_DATA%PE)         
         end if
      else 
         ierr = blas_error_param
         return
      end if
