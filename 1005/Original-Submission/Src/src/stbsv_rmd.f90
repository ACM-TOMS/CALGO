SUBROUTINE STBSV_RMD(uplo, trans, diag, n, k, A, lda, x, incx, Aa, xa, wrk, sel)
  implicit none
  character :: uplo, trans, diag
  integer :: n, k, lda, incx
  real :: A(lda,*), Aa(lda,*), x(*), xa(*), wrk(*)
  character(2) :: sel

  ! PURPOSE
  !    Calculate the reverse mode derivative of the STBSV from BLAS.
  !
  ! ARGUMENTS
  !    If STBSV was called with the following arguments:
  !
  !       uplo, trans, diag, n, k, A, lda, x, incx
  !
  !    then the corresponding call to STBSV_RMD should begin with the same arguments
  !
  !       uplo, trans, diag, n, k, A, lda, x, incx
  !
  !    with the same values. All these arguments will remain unchanged on exit.
  !    In addition the following arguments should be provided:
  !
  !    Aa
  !       (input, output, real matrix of the same dimensions as A, stored in the same
  !       band form as A, and stored in the same half according to uplo)
  !       Aa += adjoint of A due to the STBSV call
  !
  !    xa
  !       (input, output, real vector of the same dimension and increment as x)
  !       On entry: The adjoint of the x produced by STBSV
  !       On exit: The adjoint of the x supplied to STBSV
  !
  !    wrk
  !       (output, real vector of dimension at least n)
  !       When sel(2:2) = '0' so that a new xa should not be computed it is
  !       necessary to supply STBSV_RMD with a workspace vector. When sel(2:2) =
  !       '1', wrk is not referenced, because xa serves its purpose. In this case
  !       a dummy value may be given instead
  !
  !    sel
  !       (input, character*2)
  !       Used to select which adjoints to update/compute:
  !          sel(1:1) = '1' if Aa should be updated, else sel(1:1) = '0'
  !          sel(2:2) = '1' if xa should be computed, else sel(2:2) = '0'
  !       For example, if sel = '01', then only xa will be computed
  !
  ! OPERATIONS
  !    Same as for STRSV_RMD except that A and Aa use banded storage
  
  ! Local variables
  integer :: i, space, vindex, d
  logical sela, selx, ww, wx
  character :: ntrans, trno(2) = ['T', 'N']

  sela = sel(1:1) == '1'
  selx = sel(2:2) == '1'
  ntrans = trno(1 + (index('NnTt', trans)-1)/2) ! N-->T, T-->N
  wx = selx
  ww = .not.selx
  
  if(diag == 'n' .or. diag == 'N') then
    d = 0
  else
    d = 1
  end if

  ! Calculate xa (must be done before calculating Aa)
  if (selx) then
    call stbsv(uplo, ntrans, diag, n,k , A, lda, xa, incx)
  else
    call scopy(n, xa, incx, wrk, 1)
    call stbsv(uplo, ntrans, diag, n,k , A, lda, wrk, 1)
  endif

  if(sela) then
    ! Calculate Aa
    if(trans == 'n' .or. trans == 'N') then
      ! BLAS operation
      !   x = A^(-1)*x
      ! RMD operation
      !   Aa -= xa*x'
      if(uplo == 'u' .or. uplo == 'U') then
        ! A is upper triangular
        space = k
        do i = 1,n
          vindex = max(1,i-k)
          if (1+k > space+d) then
            if (wx) call saxpy(1+k-space-d, -x(1+(i-1)*incx), xa(1+(vindex-1)*incx),incx, Aa(space+1,i),1)
            if (ww) call saxpy(1+k-space-d, -x(1+(i-1)*incx), wrk(vindex)      &
              &    ,1   , Aa(space+1,i),1)
          endif
          if(space > 0) space = space - 1
        enddo
      else
        ! A is lower triangular
        space = 0

        ! Fix due to out of bounds errors emitted in this case.
        ! It will cause the index to Aa to be 2 while shape(Aa)=[1,1]. 
        ! This is not an error in this case since the number of elements 
        ! to be copied is 0.
        if(.not. (1+k-space-d==0 .and. d==1)) then
          do i = 1,n
            if (1+k > space+d) then
              if(wx) call saxpy(1+k-space-d, -x(1+(i-1)*incx), xa(1+(i+d-1)*incx),incx, Aa(1+d,i),1)
              if(ww) call saxpy(1+k-space-d, -x(1+(i-1)*incx), wrk(i+d)          &
                &,1   , Aa(1+d,i),1)
            endif
            if(i >= n-k) space = space + 1
          enddo
        endif
      endif

    else
      ! BLAS operation
      !   x = A'*x
      ! RMD operation
      !   Aa -= x*xa'
      if(uplo == 'u' .or. uplo == 'U') then
        ! A is upper triangular
        space = k
        do i = 1,n
          vindex = max(1,i-k)
          if (1+k > space+d) then
            if (wx) call saxpy(1+k-space-d, -xa(1+(i-1)*incx), x(1+(vindex-1)*incx),incx, Aa(space+1,i),1)
            if (ww) call saxpy(1+k-space-d, -wrk(i)          , x(1+(vindex-1)&
              &*incx),incx, Aa(space+1,i),1)
          endif
          if(space > 0) space = space - 1
        end do
      else
        ! A is lower triangular
        space = 0

        ! Fix due to out of bounds errors emitted in this case.
        ! See explanation above.
        if(.not. (1+k-space-d==0 .and. d==1)) then
          do i = 1,n
            if (1+k > space+d) then
              if (wx) call saxpy(1+k-space-d, -xa(1+(i-1)*incx), x(1+(i+d-1)*incx),incx, Aa(1+d,i),1)
              if (ww) call saxpy(1+k-space-d, -wrk(i)          , x(1+(i+d-1)&
                &*incx),incx, Aa(1+d,i),1)
            endif
            if(i >= n-k) space = space + 1
          end do
        end if
      end if
    end if
  end if

END SUBROUTINE STBSV_RMD
