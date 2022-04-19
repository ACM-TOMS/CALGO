SUBROUTINE SSYRK_RMD(uplo, trans, n, k, alpha, A, lda, beta, ldc, Aa, Ca, sel)
  implicit none
  character :: uplo, trans
  integer :: n, k, lda, ldc
  real :: alpha, beta
  real :: A(lda,*), Aa(lda,*), Ca(ldc,*)
  character(2) :: sel

  ! PURPOSE
  !    Calculate the reverse mode derivative of SSYRK from BLAS.
  !
  ! ARGUMENTS
  !    If SSYRK was called with the arguments
  !    
  !       uplo, trans, n, k, alpha, A, lda, beta, C, ldc
  !    
  !    then the corresponding call to SSYRK_RMD should begin with the arguments
  !    
  !       uplo, trans, n, k, alpha, A, lda, beta, ldc
  !    
  !    with the same values. All these arguments will remain unchanged on exit.
  !    Note that C is omitted. In addition the following arguments should be
  !    provided:
  !
  !    Aa
  !       (input, output, real matrix of the same dimensions as A)
  !       Aa += the adjoint of A due to the SSYRK call
  !
  !    Ca
  !       (input, output, real triangular matrix of the same dimensions as C,
  !       and stored in the same half according to uplo)
  !       On entry: the adjoint of the C produced by SSYRK
  !       On exit: the adjoint of the C supplied to SSYRK
  !
  !    sel
  !       (input, character*2)
  !       Used to select which adjoints to update/compute:
  !          sel(1:1) = '1' if Aa should be updated, else sel(1:1) = '0'
  !          sel(2:2) = '1' if Ca should be computed, else sel(2:2) = '0'
  !       For example, to update only Aa, set sel = '10'.
  !
  ! OPERATIONS
  !    (for SSYRK('L', 'N'...); C lower triangular)
  !    BLAS: C := alpha*tril(A*A') + beta*C, i.e. sym(C) := alpha*A*A' + beta*sym(C)
  !    RMD:  Aa += alpha*(Ca + Ca')*A, where Ca is value on entry
  !          Ca := beta*Ca

  ! Local variables
  integer :: i, j
  logical sela, selc

  sela = sel(1:1) == '1'
  selc = sel(2:2) == '1'
  if(sela) then
    ! Calculate Aa
    if (trans == 'n' .or. trans == 'N') then 
      ! BLAS operation
      !   C = alpha*A*A' + beta*C
      ! RMD op
      !   Aa += alpha*(Ca+diag(Ca))*A
      call ssymm('l', uplo, n, k, alpha, Ca, ldc, A, lda, 1.0, Aa, lda)
      do i = 1,n
        call saxpy(k, alpha*Ca(i,i), A(i,1), lda, Aa(i,1), lda)
      end do
    else 
      ! BLAS operation 
      !   C = alpha*A'*A + beta*C
      ! RMD operation
      !   Aa += alpha*A*(Ca+diag(Ca))
      call ssymm('r', uplo, k, n, alpha, Ca, ldc, A, lda, 1.0, Aa, lda)
      do i = 1,n
        call saxpy(k, alpha*Ca(i,i), A(1,i), 1, Aa(1,i), 1)
      end do
    end if
  end if
  
  if (selc) then
    ! Calculate Ca
    ! RMD operation
    ! Ca = beta*Ca
    do j=1,n
      if (uplo == 'u' .or. uplo == 'U') then
        call sscal(j, beta, Ca(1,j), 1)
      else
        call sscal(n-j+1, beta, Ca(j,j), 1)
      endif
    enddo
  endif
END SUBROUTINE SSYRK_RMD
