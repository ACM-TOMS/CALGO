SUBROUTINE SGBMV_RMD(trans, m, n, kl, ku, alpha, A, lda, x, incx, beta, incy,&
  & Aa, xa, ya, sel)
  real :: alpha, beta
  character :: trans
  integer :: m, n, kl, ku, lda, incx, incy
  real :: A(lda,*), x(*), Aa(lda,*), xa(*), ya(*)
  character(3) :: sel

  ! PURPOSE
  !    Calculate the reverse mode derivative of the BLAS routine SGBMV
  !
  !
  ! ARGUMENTS
  !    If SGBMV was called with the arguments
  !
  !       trans, m, n, kl, ku, alpha, A, lda, x, incx, beta, y, incy
  !
  !    then SGBMV_RMD should be called with the arguments:
  !
  !       trans, m, n, kl, ku, alpha, A, lda, x, incx, beta, incy
  !
  !    with the same values. Note that y is omitted. All these arguments will
  !    remain unchanged on exit. In addition the following arguments should
  !    be provided:
  !
  !    Aa
  !       (input, output, real matrix of the same dimensions and stored in the
  !       same band form as A)
  !       Aa += the adjoint of A due to the SGBMV call.
  !
  !    xa
  !       (input, output, real vector of the same dimension and increment as x)
  !       xa += the adjoint of x due to the SGBMV call.
  !
  !    ya
  !       (input, output, real vector of the same dimension and increment as x)
  !       ya := the adjoint of y due to the SGBMV call.
  !
  !    sel
  !       (input, character*3)
  !       Used to select which adjoints to update/compute:
  !          sel(1:1) = '1' if Aa should be updated, else sel(1:1) = '0'
  !          sel(2:2) = '1' if xa should be updated, else sel(2:2) = '0'
  !          sel(3:3) = '1' if ya should be computed, else sel(3:3) = '0'
  !       For example, to update only xa, set sel = '010'.
  !       
  ! OPERATIONS
  !    Same as for SGEMV except that A and Aa use banded storage.
  
  ! Local variables
  integer :: i, topspace, bottomspace, copysize, vpos
  logical sela, selx, sely

  sela = sel(1:1) == '1'
  selx = sel(2:2) == '1'
  sely = sel(3:3) == '1'

  topspace = ku
  bottomspace = 0
  vpos = 1

  ! If trans is 'n' or 'N' then the BLAS operation is 
  !   y = alpha*A*x + beta*y
  ! Else (if trans is 't' or 'T') the BLAS operation is 
  !   y = alpha*A'*x + beta*y

  if(sela) then
    if(trans == 'n' .or. trans == 'N') then
      !   Aa += alpha*ya*x'
      do i = 1,n
        copysize = min(kl+ku+1-topspace-bottomspace,m)
        if (copysize > 0) call saxpy(copysize, alpha*x(1+(i-1)*incx), &
          & ya(1+(vpos-1)*incy),incy, Aa(1+topspace,i),1)
        if(i >= ku+1) vpos = vpos + 1
        if(i <= ku) topspace = topspace-1
        ! Black Magic
        if(i>=max(m,n)-kl .or. (m<n .and. i>=m-kl)) bottomspace = bottomspace+1
      end do
    else
      !   Aa += alpha*x*ya'
      do i = 1,n
        copysize = min(kl+ku+1-topspace-bottomspace,m)
        if (copysize > 0) call saxpy(copysize,alpha*ya(1+(i-1)*incy), &
          & x(1+(vpos-1)*incx),incx, Aa(1+topspace,i),1)
        if(i >= ku+1) vpos = vpos + 1
        if(i <= ku) topspace = topspace-1
        ! Black Magic
        if(i>=max(m,n)-kl .or. (m<n .and. i>=m-kl)) bottomspace = bottomspace+1
      end do
    end if
  end if

  if(selx) then
    if(trans == 'n' .or. trans == 'N') then
      !   xa += alpha*A'*ya
      call sgbmv('t', m,n,kl,ku, alpha,A,lda, ya,incy, 1.0,xa,incx)
    else
      !   xa += alpha*A*ya
      call sgbmv('n', m,n,kl,ku, alpha,A,lda, ya,incy, 1.0,xa,incx)
    end if
  end if

  if(sely) then
    !   ya = beta*ya
    if(trans == 'n' .or. trans == 'N') then
      call sscal(m, beta, ya, incy)
    else
      call sscal(n, beta, ya, incy)
    end if
  end if

END SUBROUTINE SGBMV_RMD
