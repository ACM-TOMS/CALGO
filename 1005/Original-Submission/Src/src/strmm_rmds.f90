SUBROUTINE STRMM_RMDS(side, uplo, transa, diag, m, n, A, lda, B0, ldb, alphaa, Ba, wrk)
  implicit none
  character :: side, uplo, transa, diag
  integer :: m, n, lda, ldb
  real :: alphaa
  real :: A(lda,*), B0(ldb,*), Ba(ldb,*), wrk(*)
  
  ! PURPOSE
  !    Calculate the adjoint of alpha for STRMM from BLAS.
  !
  !
  ! ARGUMENTS
  !    If STRMM was called with the arguments
  !    
  !       side, uplo, transa, diag, m, n, alpha, A, lda, B, ldb
  !
  !    then the corresponding call to STRMM_RMDS should begin with the arguments
  !
  !       side, uplo, transa, diag, m, n, A, lda, B0, ldb
  !
  !    which all except B0 should have the same values as they had on the STRMM
  !    call, and B0 should have the value that B had on entry to the STRMM-call
  !    (STRMM only changes the B-argument). All these arguments except B0 will
  !    remain unchanged on exit, but B0 is used as workspace by STRMM_RMDS. Note
  !    that alpha is omitted. In addition the following arguments should be
  !    provided:
  !
  !    alphaa
  !       (input, output, real scalar)
  !       alphaa += the adjoint of alpha due to the SGEMM call.
  !
  !    Ba
  !       (input, real matrix of the same dimensions as B) 
  !       The adjoint of the B produced by SGEMM
  !
  !    wrk
  !       (output, real vector of dimension max(m,n))
  !       Workspace
  !
  ! OPERATIONS
  !    (for STRMM('L', 'L', 'N'...); A lower triangular)
  !    BLAS: B := alpha*A*B0
  !    RMD:  alphaa += vec(Ba)'*vec(A*B0)

  integer :: i
  real, external :: sdot
  character nt
  if (side == 'l' .or. side == 'L') then ! A is m by m
    do i = 1, n
      call scopy(m, B0(1,i), 1, wrk, 1)
      call strmv(uplo, transa, diag, m, A, lda, wrk, 1)
      alphaa = alphaa + sdot(m, Ba(1,i), 1, wrk, 1)
    enddo
  else
    if (transa == 'n' .or. transa == 'N') then ! A is n by n
      nt = 'T'
    else
      nt = 'N'
    endif
    do i = 1, m
      call scopy(n, B0(i,1), ldb, wrk, 1)
      call strmv(uplo, nt, diag, n, A, lda, wrk, 1)
      alphaa = alphaa + sdot(n, Ba(i,1), ldb, wrk, 1)
    enddo
  endif
END SUBROUTINE STRMM_RMDS
