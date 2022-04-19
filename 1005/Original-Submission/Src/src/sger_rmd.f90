SUBROUTINE SGER_RMD(m, n, alpha, x, incx, y, incy, lda, xa, ya, Aa, sel)
  implicit none
  real :: alpha
  integer :: m, n, incx, incy, lda
  real :: x(*), y(*), Aa(lda,*), xa(*), ya(*)
  character(2) :: sel

  ! PURPOSE
  !    Calculate the reverse mode derivative of SGER from BLAS. 
  !
  !
  ! ARGUMENTS
  !    If SGER was called with the arguments
  !
  !       m, n, alpha, x, incx, y, incy, A, lda 
  !    
  !    then the corresponding call to SGER_RMD should begin with the arguments:
  !
  !       m, n, alpha, x, incx, y, incy, lda
  !
  !    with the same values. All these arguments will remain unchanged on exit.
  !    Note that A is omitted. In addition the following arguments should be
  !    provided:
  !
  !    xa
  !       (input, output, real vector of the same dimension and increment as x)
  !       xa += the adjoint of x due to the SGER call.
  !
  !    ya
  !       (input, output, real vector of the same dimension and increment as y)
  !       ya += the adjoint of y due to the SGER call.
  !
  !    Aa
  !       (input, real matrix of the same dimensions as A)
  !       The adjoint of A.
  !
  !    sel
  !       (input, character*2)
  !       Used to select which adjoints to update:
  !          sel(1:1) = '1' if xa should be updated, else sel(1:1) = '0'
  !          sel(2:2) = '1' if ya should be updated, else sel(2:2) = '0'
  !       For example, to update only xa, set sel = '10'.
  !
  ! NOTE
  !    To compute A += alpha*x*x' one may call sger with a repeated argument,
  !    e.g. call sger(n, n, alpha, x, 1, x, 1, A, n). The adjoint of x can
  !    then be computed with:
  !       call sger_rmd(n, n, alpha, x, 1, x, 1, n, xa, dummy, Aa, '10')
  !       call sscal(n, 2.0, xa, 1).  
  !       
  ! OPERATIONS
  !    BLAS: A += alpha*x*y'
  !    RMD:  xa += alpha*Aa*y
  !          ya += alpha*Aa'*x
  !          Aa unchanged

  ! Local variables
  logical selx, sely

  selx = sel(1:1) == '1'
  sely = sel(2:2) == '1'

  if(selx) call sgemv('N', m, n, alpha, Aa, lda, y, incy, 1.0, xa, incx)
  if(sely) call sgemv('T', m, n, alpha, Aa, lda, x, incx, 1.0, ya, incy)
END SUBROUTINE SGER_RMD
