SUBROUTINE SDOT_RMD(n, x, incx, y, incy, dota, xa, ya, sel)
  implicit none
  integer :: n, incx, incy
  real :: dota, x(*), y(*), xa(*), ya(*)
  character(2) :: sel

  ! PURPOSE
  !    Calculate the reverse mode derivative of SDOT from BLAS.
  !
  ! ARGUMENTS
  !    If SDOT was called with the statement
  !
  !       dot = SDOT(n, x, incx, y, incy)
  !
  !    then SDOT_RMD should be called with the same arguments:
  !
  !       n, x, incx, y, incy
  !
  !    with the same values. These arguments will remain unchanged on exit. In
  !    addition the following arguments should be provided:
  !
  !    dota
  !       (input, real scalar)
  !       The adjoint of dot.
  !
  !    xa
  !       (input, output, real vector of the same dimension and increment as x)
  !       xa += the adjoint of x due to the SDOT call.
  !
  !    ya
  !       (input, output, real vector of the same dimension and increment as y)
  !       ya += the adjoint of y due to the SDOT call.
  !
  !    sel
  !       (input, character*2)
  !       Used to select which adjoints to update:
  !           sel(1:1) = '1' if xa should be updated, else sel(1:1) = '0'
  !           sel(2:2) = '1' if ya should be updated, else sel(2:2) = '0'
  !       For example, to update only xa, set sel = '10'.
  !
  ! NOTE
  !    To compute the adjoint of square norm, s = sdot(n, x, 1, x, 1) one may
  !    use the following calls:
  !
  !       call sdot_rmd(n, x, 1, x, 1, sa, xa, dummy, '10')
  !       call sscal(n, 2.0, xa, 1)
  !
  ! OPERATIONS
  !    BLAS: dot = x'*y
  !    RMD:  dota unchanged
  !          xa += y*dota
  !          ya += x*dota

  ! Local variables
  logical selx, sely
  selx = sel(1:1) == '1'
  sely = sel(2:2) == '1'

  if(selx) call saxpy(n, dota, y, incy, xa, incx)
  if(sely) call saxpy(n, dota, x, incx, ya, incy)
END SUBROUTINE SDOT_RMD
