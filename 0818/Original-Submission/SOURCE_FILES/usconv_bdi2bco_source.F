      call accessdata(XSP_DATA,a,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      if(XSP_DATA%FIDA=='BDI') then
         call get_infoa(XSP_DATA%INFOA ,'c',res,ierr)
         if(res.eq.COP_OF_SOURCE) then
            XSP_DATA%FIDA='BCO'
            allocate(XSP_DATA%IA2(2))
            call get_infoa(XSP_DATA%INFOA,'d',lb,ierr) !row-dim of a block
            call get_infoa(XSP_DATA%INFOA,'f',BLDA,ierr) !blocks per diagonal
            mb=floor(real(XSP_DATA%M/lb))
            kb=floor(real(XSP_DATA%K/lb))
            call set_infoa(XSP_DATA%INFOA,'f',mb,ierr) !row-dim in blocks
            call set_infoa(XSP_DATA%INFOA,'g',kb,ierr) !col-dim in blocks
            call XPRE_BDI2BCO(XSP_DATA%A, XSP_DATA%IA1,& 
                              XSP_DATA%IA2,BLDA,BNNZ,lb)
            call set_infoa(XSP_DATA%INFOA,'n',BNNZ,ierr)    
         end if
      else 
         ierr = blas_error_param
         return
      end if
