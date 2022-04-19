SUBROUTINE SPOTRF_RMD(uplo, n, A, lda, Aa)
  implicit none
  character :: uplo
  integer   :: n, lda
  real :: A(lda,*), Aa(lda,*)

  ! PURPOSE
  !    Calculate the reverse mode derivative of the Lapack Cholesky factorization
  !    subroutine SPOTRF.
  !
  ! ARGUMENTS
  !    If SPOTRF was called with the arguments
  !
  !       uplo, n, A, lda, info
  !
  !    and finished successfully, then the corresponding call to SPOTRF_RMD
  !    should begin with the arguments:
  !
  !       uplo, n, A, lda
  !
  !    with the values which they had on exit from SPOTRF. In particular A
  !    should contain the Cholesky factor of the original matrix. All these
  !    arguments will remain unchanged on exit from SPOTRF_RMD. In addition the
  !    following argument should be provided:
  !
  !    Aa
  !       (input, output, real triangular matrix of the same dimensions as A, and
  !       stored in the same half according to uplo)
  !       On entry: The adjoint of the Cholesky factor L
  !       On exit: The adjoint of the original matrix A due to the SPOTRF call.
  !
  ! OPERATIONS
  !    BLAS: A := solution to L*L' = sym(A) (i.e. A := Cholesky factor of sym(A))
  !    RMD:  Aa := adjoint of A due to the BLAS operation
  ! 
  ! NOTES
  !    1) The call to SPOTRF must have returned with info = 0
  !    2) Observe that Aa is assigned to and not added to
  !    3) On entry to SPTORF A is in the upper or the lower triangle of the
  !       parameter A and on exit the Cholesky factor L is in the same triangle.
  !       On entry to SPOTRF_RMD the parameters are:
  !          A:  Cholesky factor, L
  !          Aa: the adjoint of L
  !       and on exit:
  !          A:  unchanged
  !          Aa: the adjoint of A
  !
  ! ALGORITHM
  !    The algorithm below is obtained by finding, line by line, the adjoint
  !    of the "recursive" version of Cholesky factorization, which can be
  !    derived as follows. Consider the block matrix equalities:
  !      LL' =  | d   0  | * | d  l1' | = | d^2   d*l1'         | = | a11  a' |
  !             | l1  L1 |   | 0  L1' |   | l1*d  L1*L1'+l1*l1' |   | a    A1 |
  !
  !    From these one obtains the following formulae for d, l1 and L1:
  !      d = sqrt(a11)
  !      l1 = a/d
  !      B1 = A1 - l1*l1'
  !      L1 = chol(B1)
  !
  !    the last one of which can be applied recursively until B1 is empty.
  
  integer :: k, i, m
  real :: akk, sdot
  logical lower
  external ssymv, sscal, sdot
  lower = uplo == 'L' .or. uplo == 'l'
  do k = n, 1, -1
    m = n - k
    akk = A(k,k)
    if (lower) then
      if (m > 0) then
        call ssymv('l', m, -1.0, Aa(k+1,k+1), lda, A(k+1,k), 1, 1.0, Aa(k+1,k), 1)
        do i=k+1,n  ! subtract Dg(Ca)*v from va
          Aa(i,k) = Aa(i,k) - Aa(i,i)*A(i,k)
        enddo
        call sscal(m, 1.0/akk, Aa(k+1,k), 1)
        Aa(k,k) = Aa(k,k) - sdot(m, A(k+1,k), 1, Aa(k+1,k), 1)
      endif
    else
      if (m > 0) then
        call ssymv('u', m, -1.0, Aa(k+1,k+1), lda, A(k,k+1), lda, 1.0, Aa(k,k+1), lda)
        do i=k+1,n  ! subtract Dg(Ca)*v from va
          Aa(k,i) = Aa(k,i) - Aa(i,i)*A(k,i)
        enddo
        call sscal(m, 1.0/akk, Aa(k,k+1), lda)
        Aa(k,k) = Aa(k,k) - sdot(m, A(k,k+1), lda, Aa(k,k+1), lda)
      endif
    endif
    Aa(k,k) = Aa(k,k)/(2.0*akk)
  enddo
END SUBROUTINE SPOTRF_RMD
