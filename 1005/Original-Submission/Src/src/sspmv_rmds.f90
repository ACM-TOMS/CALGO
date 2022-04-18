SUBROUTINE SSPMV_RMDS(uplo, n, AP, x, incx, y0, incy, alphaa, betaa, ya, sel)
  implicit none
  character :: uplo
  integer :: n, incx, incy
  real :: AP(*), x(*), y0(*), ya(*), alphaa, betaa
  character(2) :: sel

  ! PURPOSE
  !    Calculate the adjoint of alpha and/or beta for SSPMV from BLAS.
  !
  !
  ! ARGUMENTS
  !    If SSPMV was called with the arguments
  !
  !       uplo, n, alpha, AP, x, incx, beta, y, incy
  !
  !    then the corresponding call to SSPMV_RMDS should begin with the same
  !    arguments
  !
  !       uplo, n, AP, x, incx, y, incy
  !
  !    which all except y0 should have the same values as they had on the SSPMV-
  !    call, and y0 should have the value that C had on entry to the SSPMV-call
  !    (SSPMV only changes the C-argument). All these arguments except C0 will
  !    remain unchanged on exit, but y0 is used as workspace by SSPMV_RMDS. Note
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
  !       The adjoint of the y produced by SSPMV
  !
  !    sel
  !       (input, character*2)
  !       Used to select which adjoints to update/compute:
  !          sel(1:1) = '1' if alphaa should be updated, else sel(1:1) = '0'
  !          sel(2:2) = '1' if betaa should be updated, else sel(2:2) = '0'
  !       For example, to update only betaa, set sel = '01'.
  !       
  ! NOTE
  !    ya must not have been updated when SSPMV_RMDS is called and therefore a
  !    potential call to SSPMV_RMD must come after a corresponding call to
  !    SSPMV_RMDS.
  !
  ! OPERATIONS
  !    (with uplo = 'L')
  !    BLAS: y = alpha*sym(AP)*x + beta*y0    AP packed lower triangular
  !    RMD:  alphaa += ya'*sym(AP)*x
  !          betaa += ya'*y0
  !    (with uplo = 'U')
  !    BLAS: y = alpha*sym(AP')*x + beta*y0   AP packed upper triangular
  !    RMD:  alphaa += ya'*sym(AP')*x
  !          betaa += ya'*y0

  real, external :: sdot

  if(sel(2:2) == '1') then
    betaa = betaa + sdot(n, ya, incy, y0, incy)
  end if

  if(sel(1:1) == '1') then
    call sspmv(uplo, n, 1.0, AP, x, incx, 0.0, y0, incy) 
    alphaa = alphaa + sdot(n, ya, incy, y0, incy)
  end if
  
END SUBROUTINE SSPMV_RMDS

