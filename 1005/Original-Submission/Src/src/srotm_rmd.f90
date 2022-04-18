SUBROUTINE SROTM_RMD(n, x, incx, y, incy, param, xa, ya, parama, sel)
  implicit none
  integer, intent(in) :: n, incx, incy
  real, intent(in) :: x(*), y(*), param(5)
  real, intent(inout) :: xa(*), ya(*), parama(5)
  character(2), intent(in) :: sel

  ! PURPOSE
  !    Calculates the reverse mode derivative of SROTM from BLAS.
  !
  ! ARGUMENTS
  !    If SROTM was called with the arguments
  !
  !       n, x, incx, y, incy, param
  !
  !    then the corresponding call to SROTM_RMD should begin with the arguments
  !
  !       n, x, incx, y, incy, param
  !
  !    with the same values as they had on exit from the SROTM-call. All these
  !    arguments will remain unchanged on exit. In addition the following
  !    arguments should be provided:
  !
  !    xa     (input, output, real vector of the same dimension and increment as x)
  !           On entry: The adjoint of the x produced by SROTM
  !           On exit: The adjoint of the x supplied to SROTM
  !
  !    ya     (input, output, real vector of the same dimension and increment as y)
  !           On entry: The adjoint of the y produced by SROTM
  !           On exit: The adjoint of the y supplied to SROTM
  !
  !    parama (input, output, real scalar)
  !           parama += the adjoint of param due to the SROTM-call. Only entries
  !           corresponding to non-fixed param entries are updated
  !
  !    sel    (input, character*3)
  !           Used to select which adjoints to update/compute:
  !             sel(1:1) = '1' if (xa,ya) should be computed, else sel(1:1) = '0'
  !             sel(2:2) = '1' if parama should be updated, else sel(2:2) = '0'
  !           For example, to update only (xa,ya), set sel = '10'.
  !
  ! OPERATIONS
  !    BLAS: [x'; y'] := H*[x'; y']
  !          where:
  !             H = [ 1  0;  0  1] if flag = -2
  !             H = [p2 p4; p3 p5] if flag = -1
  !             H = [ 1 p4; p3  1] if flag = 0
  !             H = [p2 -1;  1 p5] if flag = 1
  !             flag = param(1), pi = param(i)
  !    RMD:  [ya'; xa'] := K*[ya'; xa']
  !          Ha' += inv(H)*A
  !          where:  K = [h22 h12; h21 h22] = matrix with param = [flag p5 p3 p4 p2]
  !                  a11 = dot(x,xa)  a12 = dot(x,ya)  (xa, ya are values on entry)
  !                  a21 = dot(y,xa)  a22 = dot(y,ya)  (xa, ya are values on entry)
  !          Elements of Ha corresponding to fixed elements of H remain unchanged
  !
  ! NOTES
  !    1. The elements and structure of H is passed in param = [flag, p2, p3, p4, p5])
  !    2. H-elements that are -1, 0 or 1 are referred to as *fixed*
  !    3. For further details, see (a) Table 5 in the accompanying article [1],
  !       (b) Remark in srotm-rmd.f90, (c) The Netlib documentation of drotm and
  !       (d) The online NAG documentation of F06EQF (DROTM)
  !
  ! [1] K Jonasson et al. RMAD of BLAS Operations, ACM TOMS 2019.
 
  ! Remark: Derivation of RMD operations:
  !      [xa'; ya'] := H'*[xa'; ya'] 
  ! so that:
  !      [ya'; xa'] := K*[ya'; xa']
  ! Also:
  !      Ha += [xa'; ya']*[x y] where x, y are values supplied to SROTM
  ! or: 
  !      Ha += [xa'; ya']*[x y]*inv(H)' where x, y are values produced by SROTM
  ! or:
  !      Ha' += inv(H)*A where A = [x'; y']*[xa ya]

  real, external :: sdot
  real :: h11, h12, h21, h22, a11, a12, a21, a22, u, paramK(5)
  integer :: flag

  flag = nint(param(1))
  if (sel(2:2) == '1') then
    if (flag == -1 .or. flag == 0) then
      h21 = param(3)
      h12 = param(4)
    endif
    if (flag == -1 .or. flag == 1) then
      h11 = param(2)
      h22 = param(5)
    endif
    a11 = sdot(n, x, incx, xa, incx)
    a12 = sdot(n, x, incx, ya, incx)
    a21 = sdot(n, y, incx, xa, incx)
    a22 = sdot(n, y, incx, ya, incx)
    if (flag == -1) then
      u = h11*h22 - h12*h21
      parama(2) = parama(2) + (h22*a11 - h12*a21)/u
      parama(3) = parama(3) + (h22*a12 - h12*a22)/u
      parama(4) = parama(4) + (-h21*a11 + h11*a21)/u
      parama(5) = parama(5) + (-h21*a12 + h11*a22)/u
    elseif (flag == 0) then
      u = 1 - h12*h21
      parama(3) = parama(3) + (a12 - h12*a22)/u
      parama(4) = parama(4) + (-h21*a11 + a21)/u
    elseif (flag == 1) then
      u = h11*h22 + 1
      parama(2) = parama(2) + (h22*a11 - a21)/u
      parama(5) = parama(5) + (a12 + h11*a22)/u
    endif
  endif
  if (sel(1:1) == '1') then
    paramK(1) = param(1)
    select case(flag)
    case(-1)
      paramK(2) = param(5)
      paramK(3) = param(3)
      paramK(4) = param(4)
      paramK(5) = param(2)
    case(0)
      paramK(3) = param(3)
      paramK(4) = param(4)
    case(1)
      paramK(2) = param(5)
      paramK(5) = param(2)
    end select
    call srotm(n, ya, incy, xa, incx, paramK)
  endif
  
END SUBROUTINE SROTM_RMD
