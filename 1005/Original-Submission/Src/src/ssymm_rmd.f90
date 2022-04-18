SUBROUTINE SSYMM_RMD(side, uplo, m, n, alpha, A, lda, B, ldb, beta, ldc, Aa, Ba, Ca, sel)
  implicit none
  character :: side, uplo
  integer :: j, m, n, lda, ldb, ldc
  real :: alpha, beta
  real :: A(lda,*), B(ldb,*), Aa(lda,*), Ba(ldb,*), Ca(ldc,*)
  character(3) :: sel

  ! PURPOSE
  !    Calculate the reverse mode derivative of SSYMM from BLAS.
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
  !       side, uplo, m, n, alpha, A, lda, B, ldb, beta, ldc
  !
  !    with the same values. All these arguments will remain unchanged on exit.
  !    Note that C is omitted. In addition the following arguments should be
  !    provided:
  !
  !    Aa
  !       (input, output, real triangular matrix of the same dimensions as A, and
  !       stored in the same half according to uplo)
  !       Aa += the adjoint of A due to the SSYMM call.
  !
  !    Ba
  !       (input, output, real matrix of the same dimensions as B)
  !       Ba += the adjoint of B due to the SSYMM call.
  !
  !    Ca
  !       (input, output, real matrix of the same dimensions as C)
  !       On entry: the adjoint of the C produced by SSYMM
  !       On exit: the adjoint of the C supplied to SSYMM
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
  !    (for SSYMM('L', 'L...))
  !    BLAS: C = alpha*sym(A)*B + beta*C, where A is a lower triangular matrix
  !    RMD:  Aa += alpha*(tril(B*Ca'+Ca*B') - diag(B*Ca'))
  !          Ba += alpha*sym(A)*Ca
  !          Ca := beta*Ca
  !
  !    (for SSYMM('R', 'L...))
  !    BLAS: C = alpha*B*sym(A) + beta*C, where A is a lower triangular matrix
  !    RMD:  Aa += alpha*(tril(B'*Ca+Ca'*B) - diag(B'*Ca))
  !          Ba += alpha*sym(A)*Ca
  !          Ca := beta*Ca

  ! Local variables
  logical sela, selb, selc
  integer i
  real, external :: sdot
  sela = sel(1:1) == '1'
  selb = sel(2:2) == '1'
  selc = sel(3:3) == '1'

  if(sela) then
    if(side == 'l' .or. side == 'L') then
      ! RMD operation:
      !   Aa += tril(B*Ca' + Ca*B') - diag(B*Ca')
      ! tril-part
      call ssyr2k(uplo, 'n', m, n, alpha,B,ldb,Ca,ldc, 1.0,Aa,lda)
      ! diag-part
      do i=1,m
        Aa(i,i) = Aa(i,i) - alpha*sdot(n, B(i,1), ldb, Ca(i,1), ldc)
      enddo
    else
      ! RMD operation:
      !   Aa += tril(B'*Ca + Ca'*B) - diag(B'*Ca)
      ! tril-part
      call ssyr2k(uplo,'t', n,m, alpha,B,ldb,Ca,ldc, 1.0,Aa,lda)
      ! diag-part
      do i=1,n
        Aa(i,i) = Aa(i,i) - alpha*sdot(m, B(1,i), 1, Ca(1,i), 1)
      enddo
    end if
  end if
  
  if(selb) then
    if(side == 'l' .or. side == 'L') then
      ! RMD operation
      !   Ba += A*Ca
      call ssymm('l',uplo, m,n, alpha,A,lda,Ca,ldc, 1.0,Ba,ldb)
    else
      ! RMD operation
      !   Ba += Ca*A
      call ssymm('r',uplo, m,n, alpha,A,lda,Ca,ldc, 1.0,Ba,ldb)
    end if
  end if

  if(selc) then
    ! RMD operation:
    !   Ca = beta*Ca
    do j=1,n
      call sscal(m, beta, Ca(1,j), 1)
    enddo
    ! call sgemm('t','n',m,n,1, 0.0,0,0,0,0, beta,Ca,ldc)
    ! Let transa='t', A=B=0 and lda=ldb=0 to ensure that all memory references are valid
  end if
  
END SUBROUTINE SSYMM_RMD
