SUBROUTINE SSYR2K_RMD(uplo, trans, n, k, alpha, A, lda, B, ldb, beta, ldc, Aa, Ba, Ca, sel)
  implicit none
  character :: uplo, trans
  integer :: n, k, lda, ldb, ldc
  real :: alpha, beta
  real :: A(lda,*), B(ldb,*), Aa(lda,*), Ba(ldb,*), Ca(ldc,*)
  character(3) :: sel

  ! PURPOSE
  !    Calculate the reverse mode derivative of SSYR2K from BLAS.
  !
  ! ARGUMENTS
  !    If SSYR2K was called with the arguments
  !    
  !       uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc
  !
  !    then the corresponding call to SSYR2K_RMD should begin with the arguments
  !
  !       uplo, trans, n, k, alpha, A, lda, B, ldb, beta, ldc
  !    
  !    with the same values. All these arguments will remain unchanged on exit.
  !    Note that C is omitted. In addition the following arguments should be
  !    provided:
  !
  !    Aa
  !       (input, output, real matrix of the same dimensions as A)
  !       Aa += the adjoint of A due to the SSYR2K call
  !
  !    Ba
  !       (input, output, real matrix of the same dimensions as A)
  !       Ba += the adjoint of B due to the SSYR2K call
  !
  !    Ca
  !       (input, output, real triangular matrix of the same dimensions as C, and
  !       stored in the same half according to uplo)
  !       On entry: the adjoint of the C produced by SSYR2K
  !       On exit: the adjoint of the C supplied to SSYR2K
  !
  !    sel
  !       (input, character*3)
  !       Used to select which adjoints to update/compute:
  !          sel(1:1) = '1' if Aa should be updated, else sel(1:1) = '0'
  !          sel(2:2) = '1' if Ba should be updated, else sel(2:2) = '0'
  !          sel(3:3) = '1' if Ca should be computed, else sel(3:3) = '0'
  !       For example, to update only Aa, set sel = '100'.
  !
  ! OPERATIONS
  !    (for SSYR2K('L', 'N'...); C lower triangular)
  !    BLAS: C := alpha*tril(A*B' + B*A') + beta*C
  !    RMD:  Aa += alpha*(Ca + Ca')*B (equiv.to: Aa += alpha*(sym(Ca) + diag(Ca))*B)
  !          Ba += alpha*(Ca + Ca')*A (equiv.to: Ba += alpha*(sym(Ca) + diag(Ca))*A)
  !          Ca := beta*Ca  

  ! Local variables
  integer :: i, j
  logical sela, selb, selc

  sela = sel(1:1) == '1'
  selb = sel(2:2) == '1'
  selc = sel(3:3) == '1'
  if (sela) then
    if(trans == 'n' .or. trans == 'N') then
      call ssymm('l', uplo, n, k, alpha, Ca, ldc, B, ldb, 1.0, Aa, lda)
      do i = 1,n
        call saxpy(k, alpha*Ca(i,i), B(i,1), ldb, Aa(i,1), lda)
      end do
    else
      call ssymm('r', uplo, k, n, alpha, Ca, ldc, B, ldb, 1.0, Aa, lda)
      do i=1,n
        call saxpy(k, alpha*Ca(i,i), B(1,i), 1, Aa(1,i), 1)
      end do
    end if
  end if

  if (selb) then 
    if(trans == 'n' .or. trans == 'N') then
      call ssymm('l', uplo, n, k, alpha, Ca, ldc, A, lda, 1.0, Ba, ldb)
      do i = 1,n
        call saxpy(k, alpha*Ca(i,i), A(i,1), lda, Ba(i,1), ldb)
      end do
    else
      call ssymm('r', uplo, k, n, alpha, Ca, ldc, A, lda, 1.0, Ba, ldb)
      do i=1,n
        call saxpy(k, alpha*Ca(i,i), A(1,i), 1, Ba(1,i), 1)
      end do
    end if
  end if

  if (selc) then
    do j=1,n
      if (uplo == 'u' .or. uplo == 'U') then
        call sscal(j, beta, Ca(1,j), 1)
      else
        call sscal(n-j+1, beta, Ca(j,j), 1)
      endif
    enddo
  endif

END SUBROUTINE SSYR2K_RMD
