SUBROUTINE SDSDOT_RMD(n, x, incx, y, incy, dota, ba, xa, ya, sel)
  implicit none
  integer :: n, incx, incy
  real :: dota, ba, x(*), y(*), xa(*), ya(*)
  character(3) :: sel

  ! PURPOSE
  !    Calculate the reverse mode derivative of SDSDOT from BLAS.
  !
  ! ARGUMENTS
  !    If SDSDOT was called with the statement
  !
  !       dot = SDSDOT(n, b, x, incx, y, incy)
  !
  !    then SDSDOT_RMD should be called with the same arguments:
  !
  !       n, x, incx, y, incy
  !
  !    with the same values as they had on the SDSDOT call. These arguments will
  !    remain unchanged on exit. Note that b is omitted. In addition the
  !    following arguments should be provided:
  !
  !    dota
  !       (input, real scalar)
  !       The adjoint of dot.
  !
  !    ba
  !       (input, output, real scalar)
  !       ba += the adjoint of b due to the SDSDOT call.
  !
  !    xa
  !       (input, output, real vector of the same dimension and increment as x)
  !       xa += the adjoint of x due to the SDSDOT call.
  !
  !    ya
  !       (input, output, real vector of the same dimension and increment as y)
  !       ya += the adjoint of y due to the SDSDOT call.
  !
  !    sel
  !       (input, character*2)
  !       Used to select which adjoints to update:
  !           sel(1:1) = '1' if ba should be updated, else sel(1:1) = '0'
  !           sel(2:2) = '1' if xa should be updated, else sel(1:1) = '0'
  !           sel(3:3) = '1' if ya should be updated, else sel(2:2) = '0'
  !       For example, to update only xa, set sel = '010'.
  !
  ! NOTE
  !    SDSDOT has no double precision version and neither does SDSDOT_RMD
  !
  ! OPERATIONS
  !    BLAS: dot = b + x'*y
  !    RMD:  dota unchanged
  !          ba += dota
  !          xa += y*dota
  !          ya += x*dota

  ! Local variables
  logical selb, selx, sely
  selb = sel(1:1) == '1'
  selx = sel(2:2) == '1'
  sely = sel(3:3) == '1'
  if(selb) ba = ba + dota
  if(selx) call saxpy(n, dota, y, incy, xa, incx)
  if(sely) call saxpy(n, dota, x, incx, ya, incy)
END SUBROUTINE SDSDOT_RMD
