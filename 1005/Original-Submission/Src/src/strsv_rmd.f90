SUBROUTINE STRSV_RMD(uplo, trans, diag, n, A, lda, x, incx, Aa, xa, wrk, sel)
  implicit none
  character :: uplo, trans, diag
  integer   :: n, lda, incx
  real      :: A(lda, *), x(*), Aa(lda, *), xa(*), wrk(*)
  character :: sel*2

  ! PURPOSE
  !    Calculate the reverse mode derivative of STRSV from BLAS.
  !
  ! ARGUMENTS
  !    If STRSV was called with the arguments
  !    
  !       uplo, trans, diag, n, A, lda, x, incx
  !
  !    then the corresponding call to STRSV_RMD should begin with the same arguments
  !
  !       uplo, trans, diag, n, A, lda, x, incx
  !
  !    with the same values. All these arguments will remain unchanged on exit.
  !    In addition the following arguments should be provided:
  !
  !    Aa   (input, output, real triangular matrix of the same dimensions as A,
  !         and stored in the same half according to uplo)
  !         Aa += adjoint of A due to the STRSV call
  !
  !    xa   (input, output, real vector of the same dimension and increment as x)
  !         On entry: The adjoint of the x produced by STRSV
  !         On exit: The adjoint of the x supplied to STRSV
  !
  !    wrk  (output, real vector of dimension at least n)
  !         When sel(2:2) = '0' so that a new xa should not be computed it is
  !         necessary to supply STRSV_RMD with a workspace vector. When sel(2:2) =
  !         '1', wrk is not referenced, because xa serves its purpose. In this
  !         case a dummy value may be given instead
  !
  !    sel  (input, character*2)
  !         Used to select which adjoints to update/compute:
  !            sel(1:1) = '1' if Aa should be updated, else sel(1:1) = '0'
  !            sel(2:2) = '1' if xa should be computed, else sel(2:2) = '0'
  !         For example, if sel = '01', then only xa will be computed
  !
  ! OPERATIONS
  !    (with uplo = 'L', trans = 'N' and diag = 'N' or 'U';
  !    BLAS: x := inv(A)*x      (*), A is lower triangular
  !    RMD:  xa := inv(A)'*xa   (**)
  !          Aa -= tril(xa*x')  x is output from (*), xa is output from (**)
  !
  !    (with uplo = 'U', trans = 'N' and diag = 'N' or 'U')
  !    BLAS: x := inv(A)*x      (*), A is upper triangular
  !    RMD:  xa := inv(A)'*xa   (**)
  !          Aa -= triu(xa*x')  x is output from (*), xa is output from (**)
  !
  !    (with uplo = 'L', trans = 'T' and diag = 'N' or 'U')
  !    BLAS: x := inv(A)'*x     (*), A is lower triangular
  !    RMD:  xa := inv(A)*xa    (**)
  !          Aa -= tril(x*xa')  x is output from (*), xa is output from (**)
  !
  !    (with uplo = 'U', trans = 'T' and diag = 'N' or 'U')
  !    BLAS: x := inv(A)'*x     (*), A is upper triangular
  !    RMD:  xa := inv(A)*xa    (**)
  !          Aa -= triu(x*xa')  x is output from (*), xa is output from (**)

  character :: ntrans, trno(2) = ['T', 'N']
  
  ! Local variables
  integer :: i, d, j, jd
  logical sela, selx, upper

  ntrans = trno(1 + (index('NnTt', trans)-1)/2)
  sela = sel(1:1) == '1'
  selx = sel(2:2) == '1'
  upper = uplo== 'u' .or. uplo == 'U'

  if(diag == 'n' .or. diag == 'N') then
    d = 0
  else
    d = 1
  end if

  ! Calculate xa (must be done before calculating Aa)
  if (selx) then
    call strsv(uplo, ntrans, diag, n, A, lda, xa, incx)
  else
    call scopy(n, xa, incx, wrk, 1)
    call strsv(uplo, ntrans, diag, n, A, lda, wrk, 1)    
  endif
  if(sela) then 
    ! Calculate Aa
    if(trans == 'n' .or. trans == 'N') then
      ! BLAS operation
      !   x := A^(-1)*x where A is triangular
      ! RMD operation
      !   Aa -= tril(xa*x')
      do i = 1,n
        j = 1 + (i-1)*incx
        jd = j + d*incx
        if (upper .and. i > d) then
          if (selx) then
            call saxpy(i-d, -x(j), xa, incx, Aa(1,i), 1)
          else
            call saxpy(i-d, -x(j), wrk, 1, Aa(1,i), 1)
          endif
        endif
        if (.not.upper .and. i <= n-d) then
          if (selx) then
            call saxpy(n+1-i-d, -x(j), xa(jd), incx, Aa(i+d,i), 1)
          else
            call saxpy(n+1-i-d, -x(j), wrk(i+d), 1, Aa(i+d,i), 1)
          endif
        endif
      enddo
    else
      ! BLAS operation
      !   x := A^(-T)*x where A is triangular
      ! RMD operation
      !   Aa -= triu(x*xa')
      do i = 1,n
        j = 1 + (i-1)*incx
        jd = j + d*incx
        if (upper .and. i > d) then
          if (selx) then
            call saxpy(i-d, -xa(j), x, incx, Aa(1,i), 1)
          else
            call saxpy(i-d, -wrk(i), x, incx, Aa(1,i), 1)
          endif
        endif
        if (.not.upper .and. i <= n-d) then
          if (selx) then
            call saxpy(n+1-i-d, -xa(j), x(jd), incx, Aa(i+d,i), 1)
          else
            call saxpy(n+1-i-d, -wrk(i), x(i+d), incx, Aa(i+d,i), 1)
          endif
        endif
      enddo
    endif
  endif
END SUBROUTINE STRSV_RMD
