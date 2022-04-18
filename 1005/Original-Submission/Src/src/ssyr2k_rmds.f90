SUBROUTINE SSYR2K_RMDS(uplo, trans, n, k, A, lda, B, ldb, C0, ldc, alphaa, betaa, Ca, sel)
  implicit none
  character :: uplo, trans
  integer :: n, k, lda, ldb, ldc
  real :: A(lda,*), B(ldb,*), C0(ldc,*), Ca(ldc,*), alphaa, betaa
  character(2) :: sel

  ! PURPOSE
  !    Calculate the adjoint of alpha and/or beta for SSYR2K from BLAS
  !
  ! ARGUMENTS
  !    If SSYR2K was called with the arguments
  !    
  !       uplo, trans, n, k, alpha, A, lda, beta, C, ldc
  !    
  !    then the corresponding call to SSYR2K_RMD should begin with the arguments
  !    
  !       uplo, trans, n, k, A, lda, C0, ldc
  !    
  !    which all except C0 should have the same values as they had on the SSYR2K-
  !    call, and C0 should have the value that C had on entry to the SSYR2K-call
  !    (SSYR2K only changes the C-argument). All these arguments except C0 will
  !    remain unchanged on exit, but C0 is used as workspace by SSYR2K_RMDS. Note
  !    that alpha and beta are omitted. In addition the following arguments
  !    should be provided:
  !
  !    alphaa
  !       (input, output, real scalar)
  !       alphaa += the adjoint of alpha due to the SSYR2K call.
  !
  !    betaa
  !       (input, output, real scalar)
  !       betaa += the adjoint of beta due to the SSYR2K call.
  !
  !    Ca
  !       (input, real triangular matrix of the same dimensions as C, and stored
  !       in the same half according to uplo)
  !       The adjoint of the C produced by SSYR2K
  !
  !    sel
  !       (input, character*2)
  !       Used to select which adjoints to update/compute:
  !          sel(1:1) = '1' if alphaa should be updated, else sel(1:1) = '0'
  !          sel(2:2) = '1' if betaa should be updated, else sel(2:2) = '0'
  !       For example, to update only betaa, set sel = '01'.
  !
  ! OPERATIONS
  !    (for SSYR2K('L', 'N'...); C n by n lower triangular)
  !    BLAS: C := alpha*tril(A*B' + B*A') + beta*C
  !    RMD:  alphaa += vech(Ca)'*vech(A*B' + B*A')
  !          betaa += vech(Ca)'*vech(C)

  integer :: i
  real, external :: sdot
  
  if (sel(2:2) == '1') then ! Must use C0 before it is destroyed, if sel ='1x'
    if (uplo == 'l' .or. uplo == 'L') then
      do i = 1, n
        betaa = betaa + sdot(n+1-i, Ca(i,i), 1, C0(i,i), 1)
      enddo
    else
      do i = 1, n
        betaa = betaa + sdot(i, Ca(1,i), 1, C0(1,i), 1)
      enddo
    endif
  endif

  if(sel(1:1) == '1') then ! Update alphaa
    call ssyr2k(uplo, trans, n, k, 1.0, A, lda, B, ldb, 0.0, C0, ldc)
    if (uplo == 'l' .or. uplo == 'L') then
      do i = 1, n
        alphaa = alphaa + sdot(n+1-i, Ca(i,i), 1, C0(i,i), 1)
      enddo
    else
      do i = 1, n
        alphaa = alphaa + sdot(i, Ca(1,i), 1, C0(1,i), 1)
      enddo
    endif
  endif
END SUBROUTINE SSYR2K_RMDS
