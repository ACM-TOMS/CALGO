SUBROUTINE STRMV_RMD(uplo, trans, diag, n, A, lda, x0, incx, Aa, xa, sel)
  implicit none
  character :: uplo, trans, diag
  integer :: n, lda, incx
  real :: A(lda, *), x0(*), Aa(lda, *), xa(*)
  character(2) :: sel
  
  ! PURPOSE
  !    Calculate the reverse mode derivative of STRMV from BLAS.
  !
  !
  ! ARGUMENTS
  !    If STRMV was called with the arguments
  !
  !       uplo, trans, diag, n, A, lda, x, incx
  !
  !    then the corresponding call to STRMV_RMD should begin with the arguments
  !
  !       uplo, trans, diag, n, A, lda, x0, incx
  !
  !    which all except x0 should have the same values as they had on the STRMV
  !    call, and x0 should have the value that x had on entry to the STRMV-call
  !    (STRMV only changes the x-argument). All these arguments will remain
  !    unchanged on exit. In addition the following arguments should be
  !    provided:
  !
  !    Aa
  !       (input, output, real triangular matrix of the same dimensions as A,
  !       and stored in the same half according to uplo)
  !       Aa += the adjoint of A due to the STRMV call
  !
  !    xa
  !       (input, output, real vector of the same dimension and increment as x)
  !       On entry: the adjoint of the x produced by STRMV
  !       On exit: the adjoint of the x supplied to STRMV
  !
  !    sel
  !       (input, character*2)
  !       Used to select which adjoints to update/compute:
  !          sel(1:1) = '1' if Aa should be updated, else sel(1:1) = '0'
  !          sel(2:2) = '1' if xa should be computed, else sel(2:2) = '0'
  !       For example, if sel = '01' only xa will be computed
  !       
  ! OPERATIONS
  !    (with uplo = 'L', trans = 'N' and diag = 'N' or 'U')
  !    BLAS: x := A*x           (*), A is lower triangular
  !    RMD:  xa := A'*xa        (**)
  !          Aa += tril(xa*x')  where x is input to (*), xa is input to (**)
  !
  !    (with uplo = 'U', trans = 'N' and diag = 'N' or 'U')
  !    BLAS: x := A*x           (*), A is upper triangular
  !    RMD:  xa := A'*xa        (**)
  !          Aa += triu(xa*x')  where x is input to (*), xa is input to (**)
  !
  !    (with uplo = 'L', trans = 'T' and diag = 'N' or 'U')
  !    BLAS: x := A'*x          (*), A is lower triangular
  !    RMD:  xa := A*xa         (**)
  !          Aa += tril(x*xa')  where x is input to (*), xa is input to (**)
  !
  !    (with uplo = 'U', trans = 'T' and diag = 'N' or 'U')
  !    BLAS: x := A'*x          (*), A is upper triangular
  !    RMD:  xa := A*xa         (**)
  !          Aa += triu(x*xa')  where x is input to (*), xa is input to (**)

  ! Local variables
  integer :: i, d
  logical sela, selx
  character :: ntrans, trno(2) = ['T', 'N']

  sela = sel(1:1) == '1'
  selx = sel(2:2) == '1'
  ntrans = trno(1 + (index('NnTt', trans)-1)/2) ! N-->T, T-->N

  if(diag == 'n' .or. diag == 'N') then
    d = 0
  else
    d = 1
  end if

  if(sela) then
    ! Calculate Aa
    if(trans == 'n' .or. trans == 'N') then
      if(uplo == 'u' .or. uplo == 'U') then
        ! A is upper triangular
        do i = 1,n
          call saxpy(i-d, x0(1+(i-1)*incx), xa, incx, Aa(1,i), 1)
        end do
      else
        ! A is lower triangular
        do i = 1,n-d
          call saxpy(n-i+1-d, x0(1+(i-1)*incx), xa(1+(i+d-1)*incx), incx, &
              Aa(i+d,i), 1)
        end do
      end if
    else
      ! BLAS operation
      !   x := A'*x where A is triangular
      ! RMD operation
      !   Aa += x*xa'
      if( uplo == 'u' .or. uplo == 'U') then
        ! A is upper triangular
        do i = 1,n
          call saxpy(i-d, xa(1+(i-1)*incx), x0, incx, Aa(1,i), 1)
        end do
      else
        ! A is lower triangular
        do i = 1,n-d
          call saxpy(n-i+1-d, xa(1+(i-1)*incx), x0(1+(i+d-1)*incx), incx, &
              Aa(i+d, i), 1)
        end do
      end if
    end if
  end if

  if(selx) then
    call strmv(uplo, ntrans, diag, n, A, lda, xa, incx)
  end if
  
END SUBROUTINE STRMV_RMD

      
