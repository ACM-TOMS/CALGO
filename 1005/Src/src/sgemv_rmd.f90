SUBROUTINE SGEMV_RMD(trans, m, n, alpha, A, lda, x, incx, beta, incy, Aa, xa, ya, sel)
  real :: alpha, beta
  character :: trans
  integer :: m, n, lda, incx, incy
  real :: A(lda,*), x(*), Aa(lda,*), xa(*), ya(*)
  character(3) :: sel
  
  ! PURPOSE
  !    Calculate the reverse mode derivative of SGEMV from BLAS.
  !
  ! ARGUMENTS
  !    If SGEMV was called with the arguments
  !    
  !       trans, m, n, alpha, A, lda, x, incx, beta, y, incy
  !    
  !    then the corresponding call to SGEMV_RMD should begin with the arguments
  !
  !       trans, m, n, alpha, A, lda, x, incx, beta, incy
  !
  !    with the same values. Note that y is omitted. All these arguments will
  !    remain unchanged on exit. In addition the following arguments should
  !    be provided:
  !
  !    Aa
  !       (input, output, real matrix of the same dimensions as A and stored in
  !       the same way)
  !       Aa += the adjoint of A due to the SGEMV call.
  !
  !    xa
  !       (input, output, real vector of the same dimension and increment as x)
  !       xa += the adjoint of x due to the SGEMV call.
  !
  !    ya
  !       (input, output, real vector of the same dimension and increment as y)
  !       On entry: the adjoint of the y produced by SGEMV
  !       On exit: the adjoint of the y supplied to SGEMV
  !
  !    sel
  !       (input, character*3)
  !       Used to select which adjoints to update/compute:
  !          sel(1:1) = '1' if Aa should be updated, else sel(1:1) = '0'
  !          sel(2:2) = '1' if xa should be updated, else sel(2:2) = '0'
  !          sel(3:3) = '1' if ya should be computed, else sel(3:3) = '0'
  !       For example, to update only Aa, set sel = '100'.
  !       
  ! OPERATIONS
  !    (when trans = 'N')
  !    BLAS: y = alpha*A*x + beta*y  for general matrix A
  !    RMD:  Aa += alpha*ya*x'       where ya is value on entry
  !          xa += alpha*A'*ya       do. 
  !          ya := beta*ya
  !    (when trans = 'T')
  !    BLAS: y = alpha*A'*x + beta*y for general matrix A
  !    RMD:  Aa += alpha*x*ya'       where ya is value on entry
  !          xa += alpha*A*ya        do. 
  !          ya := beta*ya
  !
 
  ! Local variables
  logical sela, selx, sely

  sela = sel(1:1) == '1'
  selx = sel(2:2) == '1'
  sely = sel(3:3) == '1'

  ! If trans is 'n' or 'N' then the BLAS operation is
  !   y = alpha*A*x + beta*y
  ! Else (if trans is 't' or 'T') then the BLAS operation is
  !   y = alpha*A'*x + beta*y

  if(sela) then
    if(trans == 'n' .or. trans == 'N') then
      ! Aa += alpha*ya*x'
      call sger(m, n, alpha, ya, incy, x, incx, Aa, lda)
    else
      ! Aa += alpha*x*ya'
      call sger(m, n, alpha, x, incx, ya, incy, Aa, lda)
    end if
  end if

  if(selx) then
    if(trans == 'n' .or. trans == 'N') then
      ! xa += alpha*A'*ya
      call sgemv('T', m, n, alpha, A, lda, ya, incy, 1.0, xa, incx)
    else
      ! xa += alpha*A*ya
      call sgemv('N', m, n, alpha, A, lda, ya, incy, 1.0, xa, incx)
    end if
  end if

  if(sely) then
    if(trans == 'n' .or. trans == 'N') then
      ! ya = beta*ya
      call sscal(m, beta, ya, incy)
    else
      ! ya = beta*ya
      call sscal(n, beta, ya, incy)
    end if
  end if

END SUBROUTINE SGEMV_RMD
