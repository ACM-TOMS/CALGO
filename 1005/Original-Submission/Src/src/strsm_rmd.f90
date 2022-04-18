SUBROUTINE STRSM_RMD(side, uplo, transa, diag, m, n, alpha, A, lda, B, ldb, Aa, Ba, wrk, sel)
  implicit none
  character :: side, uplo, transa, diag
  integer :: m, n, lda, ldb
  real :: alpha
  real :: A(lda,*), B(ldb,*), Aa(lda,*), Ba(ldb,*), wrk(*)
  character(2) :: sel

  ! PURPOSE
  !    Calculate the reverse mode derivative of STRSM from BLAS.
  !
  !
  ! ARGUMENTS
  !    If STRSM was called with the arguments
  !
  !       side, uplo, transa, diag, m, n, alpha, A, lda, B, ldb
  !
  !    then the corresponding call to STRSM_RMD should begin with the same
  !    arguments, containing the values which they had on exit from STRSM. These
  !    arguments will remain unchanged on exit from STRSM_RMD. In addition the
  !    following arguments should be provided:
  !
  !    Aa
  !       (input, output, real matrix of the same dimensions as A)
  !       Aa += adjoint of A due to the STRSM call
  !
  !    Ba
  !       (input, output, real matrix of the same dimensions as B)
  !       On entry: The adjoint of the B produced by STRSM
  !       On exit: The adjoint of the B supplied to STRSM
  !
  !    wrk
  !       (output, real vector of dimension at least max(m,n))
  !       When sel(2:2) = '0' so that a new Ba should not be computed it is
  !       necessary to supply STRSV_RMD with a workspace vector. When sel(2:2) =
  !       '1', wrk is not referenced, because Ba serves its purpose. In this case
  !       a dummy value may be given instead
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
  !    BLAS: B := inv(A)*B      (*)
  !    RMD:  Ba := inv(A)'*Ba   (**)
  !          Aa -= tril(Ba*B')  where B and Ba are outputs from (*) and (**)

  ! Local variables
  integer :: i, d, j
  logical sela, selb
  character :: ntransa, trno(2) = ['T', 'N']

  sela = sel(1:1) == '1'
  selb = sel(2:2) == '1'
  ntransa = trno(1 + (index('NnTt', transa)-1)/2) ! N-->T, T-->N

  ! Offset for ignoring diagonal or not.
  if(diag == 'n' .or. diag == 'N') then
    d = 0
  else
    d = 1
  end if

  ! Must be done before calculating Aa
  if (.not.selb) then
    if (side=='l' .or. side=='L') then
      do j=1,n
        call strsv_rmd(uplo, transa, diag, m, A, lda, B(1,j), 1, Aa, Ba(1,j), wrk, sel)
      enddo
    else
      do j=1,m
        call strsv_rmd(uplo, ntransa, diag, n, A, lda, B(j,1), ldb, Aa, Ba(j,1), wrk, sel)
      enddo
    endif
  else
    call strsm(side, uplo, ntransa, diag, m, n, 1.0, A, lda, Ba, ldb)
    
    if(sela) then
      ! Calculate Aa
      if(side == 'l' .or. side == 'L') then
        if(transa == 'n' .or. transa == 'N') then
          ! BLAS operation
          ! B = alpha*(A^-1)*B0
          ! RMD operation
          ! Aa -= Ba*B'
          if(uplo == 'u' .or. uplo == 'U') then
            ! A is upper triangular
            do i=1,m
              call sgemv('n',i-d,n, -1.0,Ba(1,1),ldb,B(i,1),ldb, 1.0,Aa(1,i),1)
            end do
          else
            ! A is lower triangular
            do i=1,m-d
              call sgemv('n',m-i+1-d,n, -1.0,Ba(i+d,1),ldb,B(i,1),ldb, 1.0, &
                Aa(i+d,i),1)
            end do
          end if
        else
          ! BLAS operation
          ! B = alpha*(A'^-1)*B
          ! RMD operation
          ! Aa -= alpha*B*Ba'
          if(uplo == 'u' .or. uplo == 'U') then
            ! A is upper triangular
            ! NOTA sgemm...???
            do i=1,m
              call sgemv('n',i-d,n, -1.0,B(1,1),ldb,Ba(i,1),ldb, 1.0,Aa(1,i),1)
            end do
          else
            ! A is lower triangular
            do i=1,m-d
              call sgemv('n',m-i+1-d,n, -1.0,B(i+d,1),ldb,Ba(i,1),ldb, 1.0, &
                Aa(i+d,i),1)
            end do
          end if
        end if
      else
        if(transa == 'n' .or. transa == 'N') then
          ! BLAS operation
          ! B = alpha*B*(A^-1)
          ! RMD operation
          ! Aa -= B'*Ba
          if(uplo == 'u' .or. uplo == 'U') then
            ! A is upper triangular
            do i=1,n
              call sgemv('t',m,i-d, -1.0,B(1,1),ldb,Ba(1,i),1, 1.0,Aa(1,i),1)
            end do
          else
            ! A is lower triangular
            do i=1,n-d
              call sgemv('t',m,n-i+1-d, -1.0,B(1,i+d),ldb,Ba(1,i),1, 1.0,&
                Aa(i+d,i),1)
            end do
          end if
        else
          ! BLAS operation
          ! B = alpha*B*(A'^-1)
          ! RMD operation
          ! Aa -= Ba'*B
          if(uplo == 'u' .or. uplo == 'U') then
            ! A is upper triangular
            do i=1,n
              call sgemv('t',m,i-d, -1.0,Ba(1,1),ldb,B(1,i),1, 1.0,Aa(1,i),1)
            end do
          else
            ! A is lower triangular
            do i=1,n-d
              call sgemv('t',m,n-i+1-d, -1.0,Ba(1,i+d),ldb,B(1,i),1, 1.0,&
                Aa(i+d,i),1)
            end do
          end if
        end if
      end if
    end if

    ! Scale Ba with alpha
    ! call sgemm('t','n', m, n, 0, 0.0, 0.0,0, 0.0,0, alpha, Ba, ldb)
    do i=1,n
      call sscal(m, alpha, Ba(1,i), 1)
    enddo
  endif
END SUBROUTINE STRSM_RMD
