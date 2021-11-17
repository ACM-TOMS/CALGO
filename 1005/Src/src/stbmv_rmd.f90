SUBROUTINE STBMV_RMD(uplo, trans, diag, n, k, A, lda, x0, incx, Aa, xa, sel)
  implicit none
  character :: uplo, trans, diag
  integer :: n, k, lda, incx
  real :: A(lda,*), Aa(lda,*), x0(*), xa(*)
  character(2) :: sel

  ! PURPOSE
  !    Calculate the reverse mode derivative of the STBMV from BLAS.
  !
  !
  ! ARGUMENTS
  !    If STBMV was called with the following arguments:
  !
  !       uplo, trans, diag, n, k, A, lda, x, incx
  !
  !    then the corresponding call to STBMV_RMD should begin with the arguments
  !
  !       uplo, trans, diag, n, k, A, lda, x0, incx
  !
  !    which all except x0 should have the same values as they had on the STBMV
  !    call, and x0 should have the value that x had on entry to the STBMV-call
  !    (STBMV only changes the x-argument). All these arguments will remain
  !    unchanged on exit. In addition the following arguments should be
  !    provided:
  !
  !    Aa
  !       (input, output, real matrix of the same dimensions as A, stored in the same
  !       band form as A, and stored in the same half according to uplo)
  !       Aa += the adjoint of A due to the STBMV call
  !
  !    xa
  !       (input, output, real vector of the same dimension and increment as x)
  !       On entry: the adjoint of the x produced by STBMV
  !       On exit: the adjoint of the x supplied to STBMV
  !
  !    sel
  !       (input, character*2)
  !       Used to select which adjoints to update/compute:
  !          sel(1:1) = '1' if Aa should be updated, else sel(1:1) = '0'
  !          sel(2:2) = '1' if xa should be computed, else sel(2:2) = '0'
  !       For example, to compute only xa, set sel = '01'
  !       
  ! OPERATIONS
  !    Same as for STRMV_RMD except that A and Aa use banded storage
  
  ! Local variables
  integer :: i, space, vindex, d
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
      ! BLAS operation
      !   x = A*x0
      ! RMD operation
      !   Aa += xa*x0'
      if(uplo == 'u' .or. uplo == 'U') then
        ! A is upper triangular
        space = k
        do i = 1,n
          vindex = max(1,i-k)
          call saxpy(1+k-space-d, x0(1+(i-1)*incx), &
              xa(1+(vindex-1)*incx),incx, Aa(space+1,i),1)
          if(space > 0) space = space - 1
        end do
      else
        ! A is lower triangular
        space = 0
        do i = 1,n

        ! Fix due to out of bounds errors emitted in this case.
        ! It will cause the index to Aa to be 2 while shape(Aa)=[1,1]. 
        ! This is not an error in this case since the number of elements 
        ! to be copied is 0.
        if (.not. (1+k-space-d == 0 .and. d == 1)) then
          call saxpy(1+k-space-d, x0(1+(i-1)*incx), xa(1+(i+d-1)*incx),incx, &
              Aa(1+d,i),1)
          if(i >= n-k) space = space + 1
        end if

        end do
      end if

    else
      ! BLAS operation
      !   x = A'*x0
      ! RMD operation
      !   Aa += x0*xa'
      if(uplo == 'u' .or. uplo == 'U') then
        ! A is upper triangular
        space = k
        do i = 1,n
          vindex = max(1,i-k)
          call saxpy(1+k-space-d, xa(1+(i-1)*incx), &
              x0(1+(vindex-1)*incx),incx, Aa(space+1,i),1)
          if(space > 0) space = space - 1
        end do
      else
        ! A is lower triangular
        space = 0
        do i = 1,n

        ! Fix due to out of bounds errors emitted in this case.
        ! See explanation above.
        if (.not. (1+k-space-d == 0 .and. d == 1)) then
          call saxpy(1+k-space-d, xa(1+(i-1)*incx), x0(1+(i+d-1)*incx),incx, &
              Aa(1+d,i),1)
          if(i >= n-k) space = space + 1
        end if


        end do
      end if

    end if
  end if

  if(selx) then
    ! Calculate xa
    call stbmv(uplo, ntrans, diag, n,k , A, lda, xa, incx)
  end if

END SUBROUTINE STBMV_RMD
