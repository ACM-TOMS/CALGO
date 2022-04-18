SUBROUTINE SSYMM_RMDS(side, uplo, m, n, A, lda, B, ldb, C0, ldc, alphaa, betaa, Ca, sel)
  implicit none
  character :: side, uplo
  integer :: m, n, lda, ldb, ldc
  real :: A(lda,*), B(ldb,*), C0(ldc,*), Ca(ldc,*), alphaa, betaa
  character(2) :: sel

  ! PURPOSE
  !    Calculate the adjoint of alpha and/or beta for SSYMM from BLAS.
  !
  !
  ! ARGUMENTS
  !    If SSYMM was called with the arguments
  !    
  !       side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc.
  !
  !    then the corresponding call to SSYMM_RMD should begin with the same
  !    arguments
  !
  !       side, uplo, m, n, A, lda, B, ldb, C0, ldc
  !
  !    which all except C0 should have the same values as they had on the SSYMM-
  !    call, and C0 should have the value that C had on entry to the SSYMM-call
  !    (SSYMM only changes the C-argument). All these arguments except C0 will
  !    remain unchanged on exit, but C0 is used as workspace by SSYMM_RMDS. Note
  !    that alpha and beta are omitted. In addition the following arguments
  !    should be provided:
  !
  !    alphaa
  !       (input, output, real scalar)
  !       alphaa += the adjoint of alpha due to the SGEMM call.
  !
  !    betaa
  !       (input, output, real scalar)
  !       betaa += the adjoint of beta due to the SGEMM call.
  !
  !    Ca   (input, real matrix of the same dimensions as C)
  !         The adjoint of the C produced by SSYMM
  !
  !    sel
  !       (input, character*2)
  !       Used to select which adjoints to update:
  !          sel(1:1) = '1' if alphaa should be updated, else sel(1:1) = '0'
  !          sel(2:2) = '1' if betaa should be updated, else sel(2:2) = '0'
  !       For example, to update only alphaa, set sel = '10'.
  !       
  ! OPERATIONS
  !    (for SSYMM('L', 'L...))
  !    BLAS: C = alpha*sym(A)*B + beta*C, where A is a lower triangular matrix
  !    RMD:  alphaa += vech(Ca)'*vech(sym(A)*B)
  !          betaa += vech(Ca)'*vech(C0)
  !
  !    (for SSYMM('R', 'L...))
  !    BLAS: C = alpha*B*sym(A) + beta*C, where A is a lower triangular matrix
  !    RMD:  alphaa += vech(Ca)'*vech(B*sym(A))
  !          betaa += vech(Ca)'*vech(C0)

  integer i
  real, external :: sdot
  
  if(sel(2:2) == '1') then
    do i = 1, n
      betaa = betaa + sdot(m, Ca(1,i), 1, C0(1,i), 1)
    enddo
  end if

  if(sel(1:1) == '1') then
    call ssymm(side, uplo, m, n, 1.0, A, lda, B, ldb, 0.0, C0, ldc)
    do i = 1, n
      alphaa = alphaa + sdot(m, Ca(1,i), 1, C0(1,i), 1)
    enddo
  end if
  
END SUBROUTINE SSYMM_RMDS
