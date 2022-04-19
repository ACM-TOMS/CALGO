SUBROUTINE SGEMM_RMDS(transa, transb, m, n, k, A, lda, B, ldb, C0, ldc, alphaa, betaa, Ca, sel)
  character :: transa, transb
  integer :: m, n, k, lda, ldb, ldc
  real :: A(lda,*), B(ldb,*), C0(ldc,*), Ca(ldc,*), alphaa, betaa
  character(2) :: sel
  
  ! PURPOSE
  !    Calculate the adjoint of alpha and/or beta for SGEMM from BLAS.
  !
  ! ARGUMENTS
  !    If SGEMM was called with the arguments
  !    
  !       transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc
  !    
  !    then the corresponding call to SGEMM_RMDS should begin with the arguments
  !
  !       trans, transb, m, n, k, A, lda, B, ldb, C0, ldc
  !
  !    which all except C0 should have the same values as they had on the SGEMM-
  !    call, and C0 should have the value that C had on entry to the SGEMM-call
  !    (SGEMM only changes the C-argument). All these arguments except C0 will
  !    remain unchanged on exit, but C0 is used as workspace by SGEMM_RMDS. Note
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
  !         The adjoint of the C produced by SGEMM
  !
  !    sel
  !       (input, character*2)
  !       Used to select which adjoints to update:
  !          sel(1:1) = '1' if alphaa should be updated, else sel(1:1) = '0'
  !          sel(2:2) = '1' if betaa should be updated, else sel(2:2) = '0'
  !       For example, to update only alphaa, set sel = '10'.
  !
  ! NOTE
  !    Ca must not have been updated when SGEMM_RMDS is called and therefore a
  !    potential call to SGEMM_RMD must come after a corresponding call to
  !    SGEMM_RMDS.
  !
  ! OPERATIONS
  !          SGEMM('N', 'N'...)             SGEMM('T', 'N'...)
  !    BLAS: C = alpha*A*B + beta*C         C = alpha*A'*B + beta*C
  !    RMD:  alphaa += vec(Ca)'*vec(A*B)    alphaa += vec(Ca')*vec(A'*B)
  !          betaa += vec(Ca)'*vec(C0)      betaa += vec(Ca)'*vec(C0)
  !
  !          SGEMM('N', 'T'...)             SGEMM('T', 'T'...)
  !    BLAS: C = alpha*A*B' + beta*C        C = alpha*A'*B' + beta*C
  !    RMD:  alphaa += vec(Ca)'*vec(A*B')   alphaa += vec(Ca')*vec(A'*B')
  !          betaa += vec(Ca)'*vec(C0)      betaa += vec(Ca)'*vec(C0)

  integer :: i
  real, external :: sdot
    
  if(sel(2:2) == '1') then
    do i = 1, n
      betaa = betaa + sdot(m, Ca(1,i), 1, C0(1,i), 1)
    enddo
  end if

  if(sel(1:1) == '1') then
    call sgemm(transa, transb, m, n, k, 1.0, A, lda, B, ldb, 0.0, C0, ldc)
    do i = 1, n
      alphaa = alphaa + sdot(m, Ca(1,i), 1, C0(1,i), 1)
    enddo
  end if
  
END SUBROUTINE SGEMM_RMDS
