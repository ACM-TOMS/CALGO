 
    m=0
    ierr=-1
    
    if (a.ge.0) then    
       call accessdata(XSP_DATA,a,ierr)
       if (ierr.ne.0) then
          if (pname.eq.blas_valid_handle) then
             m=-1
          elseif (pname.eq.blas_invalid_handle) then
             m=1
          else
             m=-1
          end if
          return
       elseif(pname.eq.blas_valid_handle) then
          m=1
       end if
    else 
       call XACCESS_MATRIX(XMATRIX,a,ierr)
       if (ierr.ne.0) then
          if (pname.eq.blas_valid_handle) then
             m=-1
          elseif (pname.eq.blas_invalid_handle) then
             m=1
          else
             m=-1
          end if
          return
       elseif(pname.eq.blas_new_handle) then
          if (XMATRIX%new.eq.1) then
             m=1
          else
             m=0
          end if
       elseif(pname.eq.blas_open_handle) then
          if (XMATRIX%new.eq.0) then
             m=1
          else
             m=0
          end if
       else
          m=-1
          return
       end if
    end if
    if(pname.eq.blas_zero_base) then 
       call get_descra(XSP_DATA%DESCRA,'b',test,ierr)    
       if(test.eq.'C') then
          m=1
       end if
    elseif(pname.eq.blas_one_base) then 
       call get_descra(XSP_DATA%DESCRA,'b',test,ierr)    
       if(test.eq.'F') then
          m=1
       end if
    elseif(pname.eq.blas_general) then 
       call get_descra(XSP_DATA%DESCRA,'t',test,ierr)    
       if(test.eq.'G') then
          m=1
       end if
    elseif(pname.eq.blas_symmetric) then 
       call get_descra(XSP_DATA%DESCRA,'t',test,ierr)    
       if(test.eq.'S') then
          m=1
       end if
    elseif(pname.eq.blas_hermitian) then
       call get_descra(XSP_DATA%DESCRA,'t',test,ierr)    
       if(test.eq.'H') then
          m=1
       end if
    elseif(pname.eq.blas_upper_triangular) then 
       call get_descra(XSP_DATA%DESCRA,'t',test,ierr)    
       if(test.eq.'T') then
          call get_descra(XSP_DATA%DESCRA,'a',test,ierr)
          if(test.eq.'U') then
             m=1
          end if
       end if
    elseif(pname.eq.blas_lower_triangular) then
       call get_descra(XSP_DATA%DESCRA,'t',test,ierr)    
       if(test.eq.'T') then
          call get_descra(XSP_DATA%DESCRA,'a',test,ierr)
          if(test.eq.'L') then
             m=1
          end if
       end if
    elseif(pname.eq.blas_row_major) then
       call get_descra(XSP_DATA%DESCRA,'f',test,ierr)    
       if(test.eq.'R') then
          m=1
       end if
    elseif(pname.eq.blas_col_major) then
       call get_descra(XSP_DATA%DESCRA,'f',test,ierr)    
       if(test.eq.'C') then
          m=1
       end if
    elseif(pname.eq.blas_complex) then
       if ((rest.eq.CSP_MATRIX).or.(rest.eq.ZSP_MATRIX)) then
          m=1
       else
          m=0
       end if
    elseif(pname.eq.blas_real) then
       if ((rest.eq.SSP_MATRIX).or.(rest.eq.DSP_MATRIX)) then
          m=1
       else
          m=0
       end if
    elseif(pname.eq.blas_integer) then
       if (rest.eq.ISP_MATRIX) then
          m=1
       else
          m=0
       end if
    elseif(pname.eq.blas_double_precision) then
       if ((rest.eq.DSP_MATRIX).or.(rest.eq.ZSP_MATRIX)) then
          m=1
       else
          m=0
       end if      
    elseif(pname.eq.blas_single_precision) then
       if ((rest.eq.SSP_MATRIX).or.(rest.eq.CSP_MATRIX)) then
          m=1
       else
          m=0
       end if      
    elseif(pname.eq.blas_num_rows) then
       m=XSP_DATA%M   
    elseif(pname.eq.blas_num_cols) then
       m=XSP_DATA%K
    elseif(pname.eq.blas_num_nonzeros) then
       call get_infoa(XSP_DATA%INFOA,'n',m,ierr)
    else
       m=-1
       return
    end if
