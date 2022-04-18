SUBROUTINE STRSM_RMDS(side, uplo, transa, diag, m, n, A, lda, B0, ldb, alphaa, Ba)
  implicit none
  character :: side, uplo, transa, diag
  integer :: m, n, lda, ldb
  real :: alphaa
  real :: A(lda,*), B0(ldb,*), Ba(ldb,*)
  
  ! PURPOSE
  !    Calculate the adjoint of alpha for STRSM from BLAS.
  !
  !
  ! ARGUMENTS
  !    If STRSM was called with the arguments
  !    
  !       side, uplo, transa, diag, m, n, alpha, A, lda, B, ldb
  !
  !    then the corresponding call to STRSM_RMDS should begin with the arguments
  !
  !       side, uplo, transa, diag, m, n, A, lda, B0, ldb
  !
  !    which all except B0 should have the same values as they had on the STRSM
  !    call, and B0 should have the value that B had on entry to the STRSM-call
  !    (STRSM only changes the B-argument). All these arguments except B0 will
  !    remain unchanged on exit, but B0 is used as workspace by STRSM_RMDS. Note
  !    that alpha is omitted. In addition the following arguments should be
  !    provided:
  !
  !    alphaa
  !       (input, output, real scalar)
  !       alphaa += the adjoint of alpha due to the SGEMM call.
  !
  !    Ba   (input, real matrix of the same dimensions as B)
  !         The adjoint of the B produced by SGEMM
  !
  ! OPERATIONS
  !    (for STRSM('L', 'L', 'N'...); A lower triangular)
  !    BLAS: B := alpha*inv(A)*B0
  !    RMD:  alphaa += vec(Ba)'*vec(inv(A)*B0)

  integer :: i
  real, external :: sdot
  
  call strsm(side, uplo, transa, diag, m, n, 1.0, A, lda, B0, ldb)
  do i = 1, n
    alphaa = alphaa + sdot(m, Ba(1,i), 1, B0(1,i), 1)
  enddo
  
END SUBROUTINE STRSM_RMDS
