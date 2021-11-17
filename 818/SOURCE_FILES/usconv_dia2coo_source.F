      call accessdata(XSP_DATA,a,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      if(XSP_DATA%FIDA=='DIA') then
         call get_infoa(XSP_DATA%INFOA ,'c',res,ierr)
         if(res.eq.COP_OF_SOURCE) then
            XSP_DATA%FIDA='COO'
            allocate(XSP_DATA%IA2(2))
            call get_infoa(XSP_DATA%INFOA ,'d',LDA,ierr) !row-dim of val
            call get_infoa(XSP_DATA%INFOA ,'n',NNZ,ierr)
            call XPRE_DIA2COO(XSP_DATA%A, XSP_DATA%IA1, &
                 XSP_DATA%IA2,LDA,NNZ)        
         end if
      else 
         ierr = blas_error_param
         return
      end if
