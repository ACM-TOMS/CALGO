SUBROUTINE STRMM_RMD(side, uplo, transa, diag, m, n, alpha, A, lda, B0, ldb, Aa, Ba, sel)
  implicit none
  character :: side, uplo, transa, diag
  integer :: m, n, lda, ldb
  real :: alpha
  real :: A(lda,*), B0(ldb,*), Aa(lda,*), Ba(ldb,*)
  character(2) :: sel
  
  ! PURPOSE
  !    Calculate the reverse mode derivative of STRMM from BLAS.
  !
  !
  ! ARGUMENTS
  !    If STRMM was called with the arguments
  !    
  !       side, uplo, transa, diag, m, n, alpha, A, lda, B, ldb
  !
  !    then the corresponding call to STRMM_RMD should begin with the arguments
  !
  !       side, uplo, transa, diag, m, n, alpha, A, lda, B0, ldb
  !
  !    which all except B0 should have the same values as they had on the STRMM
  !    call, and B0 should have the value that B had on entry to the STRMM-call
  !    (STRMM only changes the B-argument). All these arguments will remain
  !    unchanged on exit. In addition the following arguments should be
  !    provided:
  !
  !    Aa
  !       (input, output, real matrix of the same dimensions as A) The
  !       Aa := the adjoint of A due to the STRMM call
  !
  !    Ba
  !       (input, output, real vector of the same dimension and increment as x)
  !       On entry: the adjoint of the B produced by STRMM
  !       On exit: the adjoint of the B supplied to STRMM
  !
  !    sel
  !       (input, character*2)
  !       Used to select which adjoints to update/compute:
  !          sel(1:1) = '1' if Aa should be updated, else sel(1:1) = '0'
  !          sel(2:2) = '1' if Ba should be computed, else sel(2:2) = '0'
  !       For example, to update only Aa, set sel = '10'.
  !
  ! OPERATIONS
  !    (for STRMM('L', 'L', 'N'...); A lower triangular)
  !    BLAS: B := alpha*A*B           (*)
  !    RMD:  Ba := alpha*A'*Ba        (**)
  !          Aa += alpha*tril(Ba*B')  where B and Ba are inputs to (*) and (**)
  
  ! Local variables
  integer :: i, d
  logical sela, selb
  character :: ntransa, trno(2) = ['T', 'N']

  sela = sel(1:1) == '1'
  selb = sel(2:2) == '1'
  ntransa = trno(1 + (index('NnTt', transa)-1)/2) ! N-->T, T-->N
  
  if(diag == 'n' .or. diag == 'N') then
    d = 0
  else
    d = 1
  end if

  if(sela) then
    ! Calculate Aa
    if(side == 'l' .or. side == 'L') then
      if(transa == 'n' .or. transa == 'N') then
        ! BLAS operation
        ! B = alpha*A*B0
        ! RMD operation
        ! Aa += alpha*Ba*B0'
        if(uplo == 'u' .or. uplo == 'U') then 
          ! A is upper triangular
          do i=1,m
            call sgemv('n',i-d,n, alpha,Ba(1,1),ldb,B0(i,1),ldb,1.0,Aa(1,i),1)
          end do
        else
          ! A is lower triangular
          do i=1,m-d
            call sgemv('n',m-i+1-d,n,alpha,Ba(i+d,1),ldb,B0(i,1),ldb,1.0, &
                Aa(i+d,i),1)
          end do
        end if
      else
        ! BLAS operation
        ! B = alpha*A'*B0
        ! RMD operation
        ! Aa += alpha*B*Ba'
        if(uplo == 'u' .or. uplo == 'U') then
          ! A is upper triangular
          do i=1,m
            call sgemv('n',i-d,n, alpha,B0(1,1),ldb,Ba(i,1),ldb,1.0,Aa(1,i),1)
          end do
        else
          ! A is lower triangular
          do i=1,m-d
            call sgemv('n',m-i+1-d,n,alpha,B0(i+d,1),ldb,Ba(i,1),ldb,1.0, &
                Aa(i+d,i),1)
          end do
        end if
      end if
    else
      if(transa == 'n' .or. transa == 'N') then
        ! BLAS operation
        ! B = alpha*B0*A
        ! RMD operation
        ! Aa += alpha*B0'*Ba
        if(uplo == 'u' .or. uplo == 'U') then
          ! A is upper triangular
          do i=1,n
            call sgemv('t',m,i-d,alpha,B0(1,1),ldb,Ba(1,i),1, 1.0,Aa(1,i),1)
          end do
        else
          ! A is lower triangular
          do i=1,n-d
            call sgemv('t',m,n-i+1-d,alpha,B0(1,i+d),ldb,Ba(1,i),1, 1.0, &
                Aa(i+d,i),1)
          end do
        end if
      else
        ! BLAS operation
        ! B = alpha*B0*A'
        ! RMD operation
        ! Aa += alpha*Ba'*B0
        if(uplo == 'u' .or. uplo == 'U') then
          ! A is upper triangular
          do i=1,n
            call sgemv('t',m,i-d,alpha,Ba(1,1),ldb,B0(1,i),1, 1.0,Aa(1,i),1)
          end do
        else
          ! A is lower triangular
          do i=1,n-d
            call sgemv('t',m,n-i+1-d,alpha,Ba(1,i+d),ldb,B0(1,i),1, 1.0, &
                Aa(i+d,i),1)
          end do
        end if
      end if
    end if
  end if
  
  if(selb) then
    call strmm(side,uplo,ntransa,diag,m,n,alpha,A,lda,Ba,ldb)
  end if
  
END SUBROUTINE STRMM_RMD
