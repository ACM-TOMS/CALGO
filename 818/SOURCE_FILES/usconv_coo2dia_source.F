      ierr=-1
      call accessdata(XSP_DATA,a,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      if(XSP_DATA%FIDA=='COO') then
         call get_infoa(XSP_DATA%INFOA ,'c',res,ierr)
         if(res.eq.COP_OF_SOURCE) then
            XSP_DATA%FIDA='DIA'
            call  XPRE_COO2DIA( XSP_DATA%M, XSP_DATA%K, XSP_DATA%A,& 
                 XSP_DATA%IA1, XSP_DATA%IA2,LDA,NDIAG)
            nullify(XSP_DATA%IA2)
            nnz = count( XSP_DATA%A.ne.0.)
            if(nnz.le.lda*ndiag*0.5) then
               ! Warning Many zeros stored !
            end if
            call set_infoa(XSP_DATA%INFOA,'n',nnz,ierr)
            call set_infoa(XSP_DATA%INFOA,'d',LDA,ierr) !row-dim of val
            call set_infoa(XSP_DATA%INFOA,'e',NDIAG,ierr) !col-dim of val
         end if
      else 
         ierr = blas_error_param
         return
      end if
