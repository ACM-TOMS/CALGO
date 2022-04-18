      call accessdata(XSP_DATA,a,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      if(XSP_DATA%FIDA=='BCO') then
         call get_infoa(XSP_DATA%INFOA ,'c',res,ierr)
         if(res.eq.COP_OF_SOURCE) then
            XSP_DATA%FIDA='BDI'
            call get_infoa(XSP_DATA%INFOA ,'e',lb,ierr)
            call get_infoa(XSP_DATA%INFOA ,'f',mb,ierr)
            call get_infoa(XSP_DATA%INFOA ,'g',kb,ierr)
            BLDA=min(mb,kb)
            call XPRE_BCO2BDI(mb,kb,lb,XSP_DATA%A,&
                       XSP_DATA%IA1,XSP_DATA%IA2,BLDA,BNDIAG)     
            nullify(XSP_DATA%IA2)
            call set_infoa(XSP_DATA%INFOA,'d',lb,ierr) !row-dim of a block
            call set_infoa(XSP_DATA%INFOA,'e',lb,ierr) !col-dim of a block
            call set_infoa(XSP_DATA%INFOA,'f',BLDA,ierr) !blocks per diagonal
            call set_infoa(XSP_DATA%INFOA,'g',BNDIAG,ierr) !no of diagonals
         end if
      else 
         ierr = blas_error_param
         return
      end if
