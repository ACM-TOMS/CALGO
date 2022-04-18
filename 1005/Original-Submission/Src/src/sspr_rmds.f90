SUBROUTINE SSPR_RMDS(uplo, n, x, incx, alphaa, APa)
  integer :: n, incx
  real :: alphaa, x(*), APa(*)
  character :: uplo

  ! PURPOSE
  !    Calculate the adjoint of alpha for SSPR from BLAS
  !
  !
  ! ARGUMENTS
  !    If SSPR was called with the arguments
  !
  !       uplo, n, alpha, x, incx, AP
  !    
  !    then the corresponding call to SSPR_RMD should begin with the arguments
  !
  !       uplo, n, x, incx
  !
  !    with the same values. All these arguments will remain unchanged on exit.
  !    Note that alpha and AP are omitted. In addition the following arguments
  !    should be provided:
  !
  !    alphaa
  !       (input, output, real scalar)
  !       alphaa += the adjoint of alpha due to the SSPR call.
  !
  !    APa
  !       (input, real packed triangular matrix stored in a vector with
  !       n*(n+1)/2 elements in the same way as AP)
  !       The adjoint of AP
  !
  ! OPERATIONS
  !    (with uplo = 'L')
  !    BLAS: AP += alpha*tril(x*x') i.e. sym(AP) := alpha*x*x' + sym(AP)
  !    RMD:  alphaa += x'*APa*x
  !    (with uplo = 'U')
  !    BLAS: AP += alpha*triu(x*x') i.e. sym(AP') := alpha*x*x' + sym(AP')
  !    RMD:  alphaa += x'*APa*x

  integer :: i, k, kA
  real, external :: sdot
  
  k = 1
  kA = 1
  if (uplo == 'l' .or. uplo == 'L') then
    do i = 1, n
      alphaa = alphaa + sdot(n+1-i, x(k), incx, APa(kA), 1)*x(k)
      k = k + incx
      kA = kA + n+1-i
    enddo
  else
    do i = 1, n
      alphaa = alphaa + sdot(i, x, incx, APa(kA), 1)*x(k)
      k = k + incx
      kA = kA + i
    enddo
  endif

END SUBROUTINE SSPR_RMDS
