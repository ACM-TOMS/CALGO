SUBROUTINE SSBMV_RMDS(uplo, n, k, A, lda, x, incx, y0, incy, alphaa, betaa, ya, sel)
  implicit none
  character :: uplo
  integer :: n, k, lda, incx, incy
  real :: A(lda, *), x(*), y0(*), ya(*), alphaa, betaa
  character(2) :: sel

  ! PURPOSE
  !    Calculate the adjoint of alpha and/or beta for SSBMV from BLAS.
  !
  !
  ! ARGUMENTS
  !    If SSBMV was called with the arguments
  !
  !       uplo, n, k, alpha, A, lda, x, incx, beta, y, incy
  !
  !    then the corresponding call to SSBMV_RMDS should begin with the same
  !    arguments
  !
  !       uplo, n, k, A, lda, x, incx, y, incy
  !
  !    which all except y0 should have the same values as they had on the SSBMV-
  !    call, and y0 should have the value that C had on entry to the SSBMV-call
  !    (SSBMV only changes the C-argument). All these arguments except C0 will
  !    remain unchanged on exit, but y0 is used as workspace by SSBMV_RMDS. Note
  !    that alpha and beta are omitted. In addition the following arguments
  !    should be provided:
  !
  !    alphaa
  !       (input, output, real scalar)
  !       alphaa += the adjoint of alpha due to the SGEMV call.
  !
  !    betaa
  !       (input, output, real scalar)
  !       betaa += the adjoint of beta due to the SGEMV call.
  !
  !    ya
  !       (input, real vector of the same dimension and increment as y)
  !       The adjoint of the y produced by SSBMV
  !
  !    sel
  !       (input, character*2)
  !       Used to select which adjoints to update/compute:
  !          sel(1:1) = '1' if alphaa should be updated, else sel(1:1) = '0'
  !          sel(2:2) = '1' if betaa should be updated, else sel(2:2) = '0'
  !       For example, to update only betaa, set sel = '01'.
  !       
  ! NOTE
  !    ya must not have been updated when SSBMV_RMDS is called and therefore a
  !    potential call to SSBMV_RMD must come after a corresponding call to
  !    SSBMV_RMDS.
  !
  ! OPERATIONS
  !    (with uplo = 'L')
  !    BLAS: y = alpha*sym(A)*x + beta*y0    A lower triangular band
  !    RMD:  alphaa += ya'*sym(A)*x
  !          betaa += ya'*y0
  !    (with uplo = 'U')
  !    BLAS: y = alpha*sym(A')*x + beta*y0   A upper triangular band
  !    RMD:  alphaa += ya'*sym(A')*x
  !          betaa += ya'*y0

  real, external :: sdot
  
  if(sel(2:2) == '1') then
    betaa = betaa + sdot(n, ya, incy, y0, incy)
  end if

  if(sel(1:1) == '1') then
    call ssbmv(uplo, n, k, 1.0, A, lda, x, incx, 0.0, y0, incy) 
    alphaa = alphaa + sdot(n, ya, incy, y0, incy)
  end if
  
END SUBROUTINE SSBMV_RMDS

