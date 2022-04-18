      module mod_usgp
      use mod_INSERTING
      use properties
      use representation_of_data
      contains
      subroutine usgp (a,pname,m)
      implicit none
      integer ,intent(in)::a
      integer,intent(out)::m
      integer,intent(in)::pname
      type(ispmat),pointer :: isp_data
      type(sspmat),pointer :: ssp_data
      type(dspmat),pointer :: dsp_data
      type(cspmat),pointer :: csp_data
      type(zspmat),pointer :: zsp_data
      type(i_matrix),pointer :: imatrix
      type(s_matrix),pointer :: smatrix
      type(d_matrix),pointer :: dmatrix
      type(c_matrix),pointer :: cmatrix
      type(z_matrix),pointer :: zmatrix
      integer ::rest,ierr
      character ::test
      rest = modulo(a,no_of_types)
      select case(rest)
! **********************************************************************
! **********************************************************************
      case(ISP_MATRIX)
    m=0
    ierr=-1
    if (a.ge.0) then    
       call accessdata(isp_data ,a,ierr)
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
       call iaccess_matrix (imatrix ,a,ierr)
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
          if (imatrix %new.eq.1) then
             m=1
          else
             m=0
          end if
       elseif(pname.eq.blas_open_handle) then
          if (imatrix %new.eq.0) then
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
       call get_descra(isp_data %DESCRA,'b',test,ierr)    
       if(test.eq.'C') then
          m=1
       end if
    elseif(pname.eq.blas_one_base) then 
       call get_descra(isp_data %DESCRA,'b',test,ierr)    
       if(test.eq.'F') then
          m=1
       end if
    elseif(pname.eq.blas_general) then 
       call get_descra(isp_data %DESCRA,'t',test,ierr)    
       if(test.eq.'G') then
          m=1
       end if
    elseif(pname.eq.blas_symmetric) then 
       call get_descra(isp_data %DESCRA,'t',test,ierr)    
       if(test.eq.'S') then
          m=1
       end if
    elseif(pname.eq.blas_hermitian) then
       call get_descra(isp_data %DESCRA,'t',test,ierr)    
       if(test.eq.'H') then
          m=1
       end if
    elseif(pname.eq.blas_upper_triangular) then 
       call get_descra(isp_data %DESCRA,'t',test,ierr)    
       if(test.eq.'T') then
          call get_descra(isp_data %DESCRA,'a',test,ierr)
          if(test.eq.'U') then
             m=1
          end if
       end if
    elseif(pname.eq.blas_lower_triangular) then
       call get_descra(isp_data %DESCRA,'t',test,ierr)    
       if(test.eq.'T') then
          call get_descra(isp_data %DESCRA,'a',test,ierr)
          if(test.eq.'L') then
             m=1
          end if
       end if
    elseif(pname.eq.blas_row_major) then
       call get_descra(isp_data %DESCRA,'f',test,ierr)    
       if(test.eq.'R') then
          m=1
       end if
    elseif(pname.eq.blas_col_major) then
       call get_descra(isp_data %DESCRA,'f',test,ierr)    
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
       m= isp_data %M   
    elseif(pname.eq.blas_num_cols) then
       m= isp_data %K
    elseif(pname.eq.blas_num_nonzeros) then
       call get_infoa(isp_data %INFOA,'n',m,ierr)
    else
       m=-1
       return
    end if
! **********************************************************************
! **********************************************************************
      case(SSP_MATRIX)
    m=0
    ierr=-1
    if (a.ge.0) then    
       call accessdata(ssp_data ,a,ierr)
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
       call saccess_matrix (smatrix ,a,ierr)
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
          if (smatrix %new.eq.1) then
             m=1
          else
             m=0
          end if
       elseif(pname.eq.blas_open_handle) then
          if (smatrix %new.eq.0) then
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
       call get_descra(ssp_data %DESCRA,'b',test,ierr)    
       if(test.eq.'C') then
          m=1
       end if
    elseif(pname.eq.blas_one_base) then 
       call get_descra(ssp_data %DESCRA,'b',test,ierr)    
       if(test.eq.'F') then
          m=1
       end if
    elseif(pname.eq.blas_general) then 
       call get_descra(ssp_data %DESCRA,'t',test,ierr)    
       if(test.eq.'G') then
          m=1
       end if
    elseif(pname.eq.blas_symmetric) then 
       call get_descra(ssp_data %DESCRA,'t',test,ierr)    
       if(test.eq.'S') then
          m=1
       end if
    elseif(pname.eq.blas_hermitian) then
       call get_descra(ssp_data %DESCRA,'t',test,ierr)    
       if(test.eq.'H') then
          m=1
       end if
    elseif(pname.eq.blas_upper_triangular) then 
       call get_descra(ssp_data %DESCRA,'t',test,ierr)    
       if(test.eq.'T') then
          call get_descra(ssp_data %DESCRA,'a',test,ierr)
          if(test.eq.'U') then
             m=1
          end if
       end if
    elseif(pname.eq.blas_lower_triangular) then
       call get_descra(ssp_data %DESCRA,'t',test,ierr)    
       if(test.eq.'T') then
          call get_descra(ssp_data %DESCRA,'a',test,ierr)
          if(test.eq.'L') then
             m=1
          end if
       end if
    elseif(pname.eq.blas_row_major) then
       call get_descra(ssp_data %DESCRA,'f',test,ierr)    
       if(test.eq.'R') then
          m=1
       end if
    elseif(pname.eq.blas_col_major) then
       call get_descra(ssp_data %DESCRA,'f',test,ierr)    
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
       m= ssp_data %M   
    elseif(pname.eq.blas_num_cols) then
       m= ssp_data %K
    elseif(pname.eq.blas_num_nonzeros) then
       call get_infoa(ssp_data %INFOA,'n',m,ierr)
    else
       m=-1
       return
    end if
! **********************************************************************
! **********************************************************************
      case(DSP_MATRIX)
    m=0
    ierr=-1
    if (a.ge.0) then    
       call accessdata(dsp_data ,a,ierr)
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
       call daccess_matrix (dmatrix ,a,ierr)
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
          if (dmatrix %new.eq.1) then
             m=1
          else
             m=0
          end if
       elseif(pname.eq.blas_open_handle) then
          if (dmatrix %new.eq.0) then
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
       call get_descra(dsp_data %DESCRA,'b',test,ierr)    
       if(test.eq.'C') then
          m=1
       end if
    elseif(pname.eq.blas_one_base) then 
       call get_descra(dsp_data %DESCRA,'b',test,ierr)    
       if(test.eq.'F') then
          m=1
       end if
    elseif(pname.eq.blas_general) then 
       call get_descra(dsp_data %DESCRA,'t',test,ierr)    
       if(test.eq.'G') then
          m=1
       end if
    elseif(pname.eq.blas_symmetric) then 
       call get_descra(dsp_data %DESCRA,'t',test,ierr)    
       if(test.eq.'S') then
          m=1
       end if
    elseif(pname.eq.blas_hermitian) then
       call get_descra(dsp_data %DESCRA,'t',test,ierr)    
       if(test.eq.'H') then
          m=1
       end if
    elseif(pname.eq.blas_upper_triangular) then 
       call get_descra(dsp_data %DESCRA,'t',test,ierr)    
       if(test.eq.'T') then
          call get_descra(dsp_data %DESCRA,'a',test,ierr)
          if(test.eq.'U') then
             m=1
          end if
       end if
    elseif(pname.eq.blas_lower_triangular) then
       call get_descra(dsp_data %DESCRA,'t',test,ierr)    
       if(test.eq.'T') then
          call get_descra(dsp_data %DESCRA,'a',test,ierr)
          if(test.eq.'L') then
             m=1
          end if
       end if
    elseif(pname.eq.blas_row_major) then
       call get_descra(dsp_data %DESCRA,'f',test,ierr)    
       if(test.eq.'R') then
          m=1
       end if
    elseif(pname.eq.blas_col_major) then
       call get_descra(dsp_data %DESCRA,'f',test,ierr)    
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
       m= dsp_data %M   
    elseif(pname.eq.blas_num_cols) then
       m= dsp_data %K
    elseif(pname.eq.blas_num_nonzeros) then
       call get_infoa(dsp_data %INFOA,'n',m,ierr)
    else
       m=-1
       return
    end if
! **********************************************************************
! **********************************************************************
      case(CSP_MATRIX)
    m=0
    ierr=-1
    if (a.ge.0) then    
       call accessdata(csp_data ,a,ierr)
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
       call caccess_matrix (cmatrix ,a,ierr)
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
          if (cmatrix %new.eq.1) then
             m=1
          else
             m=0
          end if
       elseif(pname.eq.blas_open_handle) then
          if (cmatrix %new.eq.0) then
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
       call get_descra(csp_data %DESCRA,'b',test,ierr)    
       if(test.eq.'C') then
          m=1
       end if
    elseif(pname.eq.blas_one_base) then 
       call get_descra(csp_data %DESCRA,'b',test,ierr)    
       if(test.eq.'F') then
          m=1
       end if
    elseif(pname.eq.blas_general) then 
       call get_descra(csp_data %DESCRA,'t',test,ierr)    
       if(test.eq.'G') then
          m=1
       end if
    elseif(pname.eq.blas_symmetric) then 
       call get_descra(csp_data %DESCRA,'t',test,ierr)    
       if(test.eq.'S') then
          m=1
       end if
    elseif(pname.eq.blas_hermitian) then
       call get_descra(csp_data %DESCRA,'t',test,ierr)    
       if(test.eq.'H') then
          m=1
       end if
    elseif(pname.eq.blas_upper_triangular) then 
       call get_descra(csp_data %DESCRA,'t',test,ierr)    
       if(test.eq.'T') then
          call get_descra(csp_data %DESCRA,'a',test,ierr)
          if(test.eq.'U') then
             m=1
          end if
       end if
    elseif(pname.eq.blas_lower_triangular) then
       call get_descra(csp_data %DESCRA,'t',test,ierr)    
       if(test.eq.'T') then
          call get_descra(csp_data %DESCRA,'a',test,ierr)
          if(test.eq.'L') then
             m=1
          end if
       end if
    elseif(pname.eq.blas_row_major) then
       call get_descra(csp_data %DESCRA,'f',test,ierr)    
       if(test.eq.'R') then
          m=1
       end if
    elseif(pname.eq.blas_col_major) then
       call get_descra(csp_data %DESCRA,'f',test,ierr)    
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
       m= csp_data %M   
    elseif(pname.eq.blas_num_cols) then
       m= csp_data %K
    elseif(pname.eq.blas_num_nonzeros) then
       call get_infoa(csp_data %INFOA,'n',m,ierr)
    else
       m=-1
       return
    end if
! **********************************************************************
! **********************************************************************
      case(ZSP_MATRIX)
    m=0
    ierr=-1
    if (a.ge.0) then    
       call accessdata(zsp_data ,a,ierr)
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
       call zaccess_matrix (zmatrix ,a,ierr)
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
          if (zmatrix %new.eq.1) then
             m=1
          else
             m=0
          end if
       elseif(pname.eq.blas_open_handle) then
          if (zmatrix %new.eq.0) then
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
       call get_descra(zsp_data %DESCRA,'b',test,ierr)    
       if(test.eq.'C') then
          m=1
       end if
    elseif(pname.eq.blas_one_base) then 
       call get_descra(zsp_data %DESCRA,'b',test,ierr)    
       if(test.eq.'F') then
          m=1
       end if
    elseif(pname.eq.blas_general) then 
       call get_descra(zsp_data %DESCRA,'t',test,ierr)    
       if(test.eq.'G') then
          m=1
       end if
    elseif(pname.eq.blas_symmetric) then 
       call get_descra(zsp_data %DESCRA,'t',test,ierr)    
       if(test.eq.'S') then
          m=1
       end if
    elseif(pname.eq.blas_hermitian) then
       call get_descra(zsp_data %DESCRA,'t',test,ierr)    
       if(test.eq.'H') then
          m=1
       end if
    elseif(pname.eq.blas_upper_triangular) then 
       call get_descra(zsp_data %DESCRA,'t',test,ierr)    
       if(test.eq.'T') then
          call get_descra(zsp_data %DESCRA,'a',test,ierr)
          if(test.eq.'U') then
             m=1
          end if
       end if
    elseif(pname.eq.blas_lower_triangular) then
       call get_descra(zsp_data %DESCRA,'t',test,ierr)    
       if(test.eq.'T') then
          call get_descra(zsp_data %DESCRA,'a',test,ierr)
          if(test.eq.'L') then
             m=1
          end if
       end if
    elseif(pname.eq.blas_row_major) then
       call get_descra(zsp_data %DESCRA,'f',test,ierr)    
       if(test.eq.'R') then
          m=1
       end if
    elseif(pname.eq.blas_col_major) then
       call get_descra(zsp_data %DESCRA,'f',test,ierr)    
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
       m= zsp_data %M   
    elseif(pname.eq.blas_num_cols) then
       m= zsp_data %K
    elseif(pname.eq.blas_num_nonzeros) then
       call get_infoa(zsp_data %INFOA,'n',m,ierr)
    else
       m=-1
       return
    end if
! **********************************************************************
! **********************************************************************
      case default
         return
      end select
      end subroutine usgp
      end module mod_usgp
