SUBROUTINE SGEMV_RMDS(trans, m, n, A, lda, x, incx, y0, incy, alphaa, betaa, ya, sel)
  character :: trans
  integer :: m, n, lda, incx, incy
  real :: A(lda,*), x(*), y0(*), ya(*), alphaa, betaa
  character(2) :: sel
  
  ! PURPOSE
  !    Calculate the adjoint of alpha and/or beta for SGEMV from BLAS.
  !
  ! ARGUMENTS
  !    If SGEMV was called with the arguments
  !    
  !       trans, m, n, alpha, A, lda, x, incx, beta, y, incy
  !    
  !    then the corresponding call to SGEMV_RMDS should begin with the arguments
  !
  !       trans, m, n, A, lda, x, incx, y0, incy
  !
  !    which all except y0 should have the same values as they had on the SGEMV-
  !    call, and y0 should have the value that y had on entry to the SGEMV-call
  !    (SGEMV only changes the y-argument). All these arguments except y0 will
  !    remain unchanged on exit, but y0 is used as workspace by SGEMV_RMDS. Note
  !    that alpha and beta are omitted. In addition the following arguments
  !    should be provided:
  !
  !    alphaa
  !       (input, output, real scalar)
  !       alphaa += the adjoint of alpha due to the SGEMV-call.
  !
  !    betaa
  !       (input, output, real scalar)
  !       betaa += the adjoint of beta due to the SGEMV-call.
  !
  !    ya
  !       (input, real vector of the same dimension and increment as y)
  !       The adjoint of the y produced by SGEMV
  !
  !    sel
  !       (input, character*2)
  !       Used to select which adjoints to update:
  !          sel(1:1) = '1' if alphaa should be updated, else sel(1:1) = '0'
  !          sel(2:2) = '1' if betaa should be updated, else sel(2:2) = '0'
  !       For example, to update only alphaa, set sel = '10'.
  !       
  ! NOTE
  !    ya must not have been updated when SGEMV_RMDS is called and therefore a
  !    potential call to SGEMV_RMD must come after a corresponding call to
  !    SGEMV_RMDS.
  !
  ! OPERATIONS
  !    (when trans = 'N')
  !    BLAS: y = alpha*A*x + beta*y0
  !    RMD:  alphaa += ya'*A*x
  !          betaa += ya'*y0
  !    (when trans = 'T')
  !    BLAS: y = alpha*A'*x + beta*y0
  !    RMD:  alphaa += ya'*A'*x
  !          betaa += ya'*y0

  real, external :: sdot
  
  if(sel(2:2) == '1') then
    if(trans == 'n' .or. trans == 'N') then
      betaa = betaa + sdot(m, ya, incy, y0, incy)
    else
      betaa = betaa + sdot(n, ya, incy, y0, incy)
    end if
  end if

  if(sel(1:1) == '1') then
    call sgemv(trans, m, n, 1.0, A, lda, x, incx, 0.0, y0, incy) 
    if(trans == 'n' .or. trans == 'N') then
      alphaa = alphaa + sdot(m, ya, incy, y0, incy)
    else
      alphaa = alphaa + sdot(n, ya, incy, y0, incy)
    end if
  end if
  
END SUBROUTINE SGEMV_RMDS
