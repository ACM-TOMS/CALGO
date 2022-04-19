SUBROUTINE SROTMG_RMD(d1, d2, x1, param, d1a, d2a, x1a, y1a, parama)
  implicit none
  real :: d1, d2, x1, param(5), d1a, d2a, x1a, y1a, parama(*)

  ! PURPOSE
  !    Calculates the reverse mode derivative of SROTMG from BLAS.
  !
  ! ARGUMENTS
  !    If SROTMG was called with the arguments
  !
  !       d1, d2, x1, y1, param
  !
  !    then the corresponding call to SROTMG_RMD should begin with arguments
  !
  !       d1, d2, x1, param
  !
  !    having the same values as they had on exit from SROTMG. These arguments
  !    will remain unchanged on exit from SROTMG_RMD. Note that y1 is omitted.
  !    In addition the following arguments should be provided:
  !
  !       d1a, d2a
  !          (input, output, real scalars)
  !          On entry: The adjoints of the d1 and d2 produced by SROTMG
  !          On exit: The adjoints of the d1 and d2 supplied to SROTMG
  !
  !       x1a
  !          (input, output, real scalar)
  !          On entry: The adjoint of the x1 produced by SROTMG
  !          On exit: The adjoint of the x1 supplied to SROTMG
  !
  !       y1a
  !          (output, real scalar)
  !          The adjoint of the y1 supplied to SROTMG
  !
  !       parama
  !          (input, output, real vector of dimension 5 or 8)
  !          On entry:
  !               The entries corresponding to elements in param with elements
  !               of the matrix H should contain the adjoints of these H
  !               elements. If parama(1) = 2 then parama should have dimension 8
  !          On exit with parama(1) = 2:
  !               Information about the computations, cf. [1]:
  !               parama(6): The value of Flag before scaling
  !               parama(7): Gamma for d1
  !               parama(8): Gamma for d2
  !
  ! NOTES
  !    A sel parameter is not offered. The adjoints of d1, d2, x1 and y1 are
  !    always computed together.
  !
  ! OPERATIONS
  !    BLAS: Provided by the formulae in [1]
  !    RMD:  Obtained by differentiating the formulae in [1]
  !
  !    See also comments in srotm_rmd.f90
  !
  ! [1] CL Lawson et. al., Basic linear algebra subprograms for Fortran usage,
  !     ACM TOMS 5, 1979, 308-323.

  real :: h11, h12, h21, h22, h11a, h12a, h21a, h22a, u, x1r, y1r, d1r, d2r
  real :: gam1, gam2, d1u, d2u, ua, p1, p2, p1a, p2a, temp, tempa, rflag
  integer :: flag
  flag = nint(param(1))

  if (parama(1) == 2.0) parama(6:8) = [param(1), 1.0, 1.0]

  if (flag == -2) then
    y1a = 0;
    return
  endif
  if (flag == -1 .or. flag == 0) then
    h21 = param(3)
    h12 = param(4)
    h21a = parama(3)
    h12a = parama(4)
  endif
  if (flag == -1 .or. flag == 1) then
    h11 = param(2)
    h22 = param(5)
    h11a = parama(2)
    h22a = parama(5)
  endif

  d1u = d1
  d2u = d2
  if (flag == -1) then
    ! Use D = H'*D~*H and H*[x1;y1] = [x1~;0] (where D~ and x1~ are values on
    ! exit from srotmg) to reconstruct D, x1 and y1 as they were on entry to srotmg
    ! Also reconstruct the scales, gam1 and gam2, applied to d1 and d2 by srotmg
    ! if such scaling was done.
    u = h11*h22 - h12*h21
    x1r = h22*x1/u
    y1r = -h21*x1/u
    d1r = h11**2*d1 + h21**2*d2
    d2r = h12**2*d1 + h22**2*d2
    ! Unscale:
    if (abs(d1r*x1r**2) > abs(d2r*y1r**2)) then
      flag = 0
      gam1 = 1/h11
      gam2 = 1/h22
      h12 = h12*gam1
      h21 = h21*gam2
      h12a = h12a/gam1
      h21a = h21a/gam2
      u = 1 - h12*h21
    else
      flag = 1
      gam1 = 1/h12
      gam2 = -1/h21
      h11 = h11*gam1
      h22 = h22*gam2
      h11a = h11a/gam1
      h22a = h22a/gam2
      u = 1 + h11*h22
    endif
    x1a = x1a/gam1
    d1u = d1/gam1**2
    d2u = d2/gam2**2
    d1a = d1a*gam1**2
    d2a = d2a*gam2**2
    rflag = flag
    if (parama(1) == 2.0) parama(6:8) = [rflag, gam1, gam2]
  elseif (flag == 0) then
    ! Reconstruct d1, d2, x1 and y1 as they were on entry to srotmg
    u = 1 - h12*h21
    d1r = d1*u
    d2r = d2*u
    x1r = x1/u
    y1r = -h21*x1r
  else ! flag == 1
    ! Reconstruct d1, d2, x1 and y1 as they were on entry to srotmg
    u = h11*h22 + 1
    d1r = d2*u
    d2r = d1*u
    y1r = x1/u
    x1r = h22*y1r
  endif

  ! Computes adjoint of pre-scaling part of srotmg
  if (flag == 0) then
    p1 = d1r*x1r
    ua = x1r*x1a
    x1a = u*x1a
    d2a = d2a/u
    ua = ua - d2u*d2a
    d1a = d1a/u
    ua = ua - d1u*d1a
    h12a = h12a - h21*ua
    h21a = h21a - h12*ua
    y1a = -h21a/x1r
    x1a = x1a + h21*y1a
    p2a = h12a/p1
    p1a = -h12*p2a
  elseif (flag == 1) then
    p2 = d2r*y1r
    y1a = u*x1a
    ua = y1r*x1a
    temp = d1u
    tempa = d1a
    d1a = d2a/u
    ua = ua - d2u*d1a
    d2a = tempa/u
    ua = ua - temp*d2a
    h11a = h11a + h22*ua
    h22a = h22a + h11*ua
    x1a = h22a/y1r
    y1a = y1a - h22*x1a
    p1a = h11a/p2
    p2a = -h11*p1a
  else
    stop 'Illegal value of flag'
  endif
  d1a = d1a + x1r*p1a
  d2a = d2a + y1r*p2a
  x1a = x1a + d1r*p1a
  y1a = y1a + d2r*p2a
END SUBROUTINE SROTMG_RMD
