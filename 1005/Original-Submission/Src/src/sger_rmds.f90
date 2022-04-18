SUBROUTINE SGER_RMDS(m, n, x, incx, y, incy, lda, alphaa, Aa)
  implicit none
  integer :: m, n, incx, incy, lda
  real :: x(*), y(*), Aa(lda,*), alphaa

  ! PURPOSE
  !    Calculate the adjoint of alpha for SGER from BLAS. 
  !
  !
  ! ARGUMENTS
  !    If SGER was called with the arguments
  !
  !       m, n, alpha, x, incx, y, incy, A, lda 
  !    
  !    then the corresponding call to SGER_RMD should begin with the arguments:
  !
  !       m, n, x, incx, y, incy, lda
  !
  !    with the same values. All these arguments will remain unchanged on exit.
  !    Note that alpha and A are omitted. In addition the following arguments
  !    should be provided:
  !
  !    alphaa
  !       (input, output, real scalar)
  !       alphaa += adjoint of alpha due to the SGER-call
  !
  !    Aa
  !       (input, real matrix of the same dimensions as A)
  !       The adjoint of A.
  !
  ! OPERATIONS
  !    BLAS: A += alpha*x*y'
  !    RMD:  alphaa += x'*Aa*y

  integer :: i, k
  real, external :: sdot
  k = 1
  do i = 1, n
    alphaa = alphaa + sdot(m, x, incx, Aa(1,i), 1)*y(k)
    k = k + incy
  enddo

END SUBROUTINE SGER_RMDS
