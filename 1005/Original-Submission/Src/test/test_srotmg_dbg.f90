PROGRAM TEST_SROTMG
  integer :: i
  real :: d0(2), H(2,2), param(5), x(2), y(2), G(2,2), d(2), x0(2) = [3, 2]
  !
  ! There is a bug in the reference BLAS srotmg and drotmg at least LAPACK
  ! versions 3.2 (Nov. 2008) until 3.8.0. The bug is not present in the
  ! original level 1 BLAS. The following lines will detect the bug:
  !
  d = [4.0, 9.0]
  x = x0
  print *,'d,x=',d,x
  call srotmg(d(1), d(2), x(1), x(2), param)
  print *,'(1)'
  call par2h(param, H)
  y = matmul(H, x0)
  call assert(approx_eqs(y(2), 0.0), 'y(2) not 0 when d = [4,9]')
  !
  d = [1e-8, 1e-8]
  x = x0
  print *,'d,x=',d,x
  call srotmg(d(1), d(2), x(1), x(2), param)
  print *,'(2)'
  call par2h(param, H)
  y = matmul(H, x0)
  call assert(approx_eqs(y(2), 0.0), 'y(2) not 0 when D = 1e-8*I')
  !
  d = [1e-24, 1e-24]
  x = x0
  print *,'d,x=',d,x
  call srotmg(d(1), d(2), x(1), x(2), param)
  print *,'(3)'
  call par2h(param, H)
  y = matmul(H, x0)
  call assert(approx_eqs(y(2), 0.0), 'y(2) not 0 when D = 1e-24*I')
  !
  do i = 1, 4
    print *,'i=',i
    call random_number(x)
    d0 = d
    x0 = x
    call srotmg(d(1), d(2), x(1), x(2), param)
    call par2h(param, H)
    y = matmul(H, x0)
    call assert(approx_eqs(y(1), x(1)), 'Incorrect x(1)')
    call assert(approx_eqs(y(2), 0.0), 'y(2) not 0')
    ! Compute Givens matrix corresponding to H
    G = H
    G(1, 1:2) = H(1,1:2)/sqrt(d0(1))
    G(2, 1:2) = H(2,1:2)/sqrt(d0(2))
    G(1:2, 1) = sqrt(d(1))*G(1:2, 1)
    G(1:2, 2) = sqrt(d(2))*G(1:2, 2)
    call assert(isorthogonal(G), 'G is not orthogonal')
  enddo
    
CONTAINS
  logical function isorthogonal(A)
    ! Return true if A is (approximately) orthogonal
    real :: A(:,:)
    real, allocatable :: ATA(:,:), I(:,:)
    integer :: n, j
    ATA = matmul(transpose(A), A)
    n = size(A,1)
    allocate(I(n,n))
    I = 0
    forall (j=1:n) I(j,j) = 1
    isorthogonal = approx_eqm(ATA, I)
  end function isorthogonal

  subroutine assert(p, msg)
    logical :: p
    character(*) :: msg
    if (.not. p) then
      print '("Error detected by test_srotmg")'
      print '(A)', msg
      stop 1
    endif
  end subroutine assert

  subroutine par2h(p, H)
    real :: p(5), H(2,2)
    select case(nint(p(1)))
    case (1)
      H = reshape([p(2), -1.0, 1.0, p(5)], [2,2])
    case(0)
      H = reshape([1.0, p(3), p(4), 1.0], [2,2])
    case(-1)
      H = reshape([p(2), p(3), p(4), p(5)], [2,2])
    case(-2)
      H = reshape([1.0, 0.0, 0.0, 1.0], [2,2])
    end select
  end subroutine par2h

  logical function approx_eqs(x, y) result(aeq)
    real :: x, y, maxv, tol = epsilon(0.0)*1e2
    maxv = max(1.0, abs(x), abs(y))
    aeq = abs(x - y) < maxv*tol
  end function approx_eqs
    
  logical function approx_eqm(x, y) result(aeq)
    real :: x(:,:), y(:,:), maxv, tol = epsilon(0.0)*1e2
    maxv = max(1.0, maxval(abs(x)), maxval(abs(y)))
    aeq = all(abs(x - y) < maxv*tol)
  end function approx_eqm

END PROGRAM TEST_SROTMG
