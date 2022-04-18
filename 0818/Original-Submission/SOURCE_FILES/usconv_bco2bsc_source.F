      call accessdata(XSP_DATA,a,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      if(XSP_DATA%FIDA=='BCO') then
         call get_infoa(XSP_DATA%INFOA ,'c',res,ierr)
         if(res.eq.COP_OF_SOURCE) then
            XSP_DATA%FIDA='BSC'
            call get_infoa(XSP_DATA%INFOA,'g',col_dim_in_blocks,ierr)
            call get_infoa(XSP_DATA%INFOA,'e',col_dim_of_block,ierr)
            allocate(XSP_DATA%PB(col_dim_in_blocks))
            allocate(XSP_DATA%PE(col_dim_in_blocks))
            call XPRE_BCO2BSC(XSP_DATA%A,XSP_DATA%IA1,&
                       XSP_DATA%IA2,col_dim_in_blocks,&
                       col_dim_of_block,XSP_DATA%PB,XSP_DATA%PE)
            nullify(XSP_DATA%IA2)
         end if
      else 
         ierr = blas_error_param
         return
      end if
