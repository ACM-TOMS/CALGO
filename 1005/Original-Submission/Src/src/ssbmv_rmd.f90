SUBROUTINE SSBMV_RMD(uplo, n, k, alpha, A, lda, x, incx, beta, incy, Aa, xa, ya, sel)
  implicit none
  character :: uplo
  integer :: n, k, lda, incx, incy
  real :: alpha, beta
  real :: A(lda, *), x(*), Aa(lda, *), xa(*), ya(*)
  character(3) :: sel

  ! PURPOSE
  !    Calculate the reverse mode derivative of SSBMV from BLAS.
  !
  !
  ! ARGUMENTS
  !    If SSBMV was called with the arguments
  !
  !       uplo, n, k, alpha, A, lda, x, incx, beta, y, incy
  !
  !    then the corresponding SSBMV_RMD call should begin with the arguments
  !
  !       uplo, n, k, alpha, A, lda, x, incx, beta, incy
  !
  !    with the same values. All these arguments will remain unchanged on exit.
  !    Note that y is omitted. In addition the following arguments should be
  !    provided:
  !
  !    Aa
  !       (input, output, real triangular band matrix of the same dimensions as
  !       A, stored in the same band form, and stored in the same half according
  !       to uplo)
  !       Aa += adjoint of A due to the SSBMV call
  !
  !    xa
  !       (input, output, real vector of the same dimension and increment as x)
  !       xa += adjoint of x due to the SSBMV call
  !
  !    ya
  !       (input, output, real vector of the same dimension and increment as y)
  !       On entry: the adjoint of the y produced by SSBMV
  !       On exit: the adjoint of the y supplied to SSBMV
  !
  !    sel
  !       (input, character*3)
  !       Used to select which adjoints to update/compute:
  !          sel(1:1) = '1' if Aa should be updated, else sel(1:1) = '0'
  !          sel(2:2) = '1' if xa should be updated, else sel(2:2) = '0'
  !          sel(3:3) = '1' if ya should be computed, else sel(3:3) = '0'
  !       For example, to update only Aa, set sel = '100'.
  !       
  ! OPERATIONS
  !    Same as for SSYMV_RMD except that A and Aa use banded storage
  
  ! Local variables
  integer :: i, space, vindex
  logical :: sela, selx, sely

  sela = sel(1:1) == '1'
  selx = sel(2:2) == '1'
  sely = sel(3:3) == '1'
  
  if(sela) then
    ! Update Aa
    if(uplo == 'u' .or. uplo == 'U') then
      ! A is upper triangular
      space = k
      do i = 1,n
        vindex = max(1, i-k)
        call saxpy(1+k-space, alpha*ya(1+(i-1)*incy), &
            x(1+(vindex-1)*incx),incx, Aa(space+1,i),1)
        call saxpy(1+k-space, alpha*x(1+(i-1)*incx), &
            ya(1+(vindex-1)*incy),incy, Aa(space+1,i),1)
        if(space > 0) space = space - 1
        Aa(k+1,i) = Aa(k+1,i) - alpha*x(1 + (i-1)*incx)*ya(1 + (i-1)*incy)
      enddo
    else
      ! A is lower triangular
      space = 0
      do i = 1,n
        call saxpy(1+k-space, alpha*ya(1+(i-1)*incy), x(1+(i-1)*incx),incx, &
            Aa(1,i),1)
        call saxpy(1+k-space, alpha*x(1+(i-1)*incx), ya(1+(i-1)*incy),incy, &
            Aa(1,i),1)
        if(i >= n-k) space = space + 1
        Aa(1,i) = Aa(1,i) - alpha*x(1 + (i-1)*incx)*ya(1 + (i-1)*incy)
      end do
    end if
  end if
  
  if(selx) then
    ! Update xa
    call ssbmv(uplo, n, k, alpha, A, lda, ya, incy, 1.0, xa, incx)
  end if

  if(sely) then
    ! Update ya
    call sscal(n, beta, ya, incy)
  end if
END SUBROUTINE SSBMV_RMD

