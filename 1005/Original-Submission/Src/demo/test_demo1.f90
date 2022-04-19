SUBROUTINE F(testcase, u, v)
  ! A test function for demo1_sumsq, check the adjoint.
  use rmd_dtesttools
  !use dispmodule
  implicit none
  type(fdata), intent(inout) :: testcase
  double precision, intent(in) :: u(*)
  double precision, intent(out) :: v(*)
  double precision, allocatable :: k(:), f0(:), FF(:,:), A0(:,:), A(:,:,:),&
    & C(:,:), d(:), ka(:), w(:)
  integer m, n, s, n1, n2, n3, n4, n5, n6
  double precision :: theta
  m = testcase % integers(1)
  n = testcase % integers(2)
  s = testcase % integers(3)
  n1 = n
  n2 = n1 + n*m
  n3 = n2 + n*n;
  n4 = n3 + n*n*m
  n5 = n4 + s*n
  n6 = n5 + s
  allocate(w(n6), k(m), f0(n), FF(n,m), A0(n,n), A(n,n,m), C(s,n), d(s), ka(m))
  w = testcase % reals
  k = u(1:m)
  f0 =         w(   1:n1)
  FF = reshape(w(n1+1:n2), [n,m])
  A0 = reshape(w(n2+1:n3), [n,n])
  A  = reshape(w(n3+1:n4), [n,n,m])
  C  = reshape(w(n4+1:n5), [s,n])
  d  =         w(n5+1:n6)
  call demo1_sumsq(m, n, s, k, f0, FF, A0, A, C, d, theta, ka)
  v(1) = theta
END SUBROUTINE F

SUBROUTINE F_RMD(testcase, u, v, ua, va)
  use rmd_dtesttools
  implicit none
  type(fdata), intent(inout) :: testcase
  double precision, intent(in) :: u(*), v(*), va(*)
  double precision, intent(inout) :: ua(*)
  double precision, allocatable :: k(:), f0(:), F(:,:), A0(:,:), A(:,:,:), C(:,:),&
    &                              d(:), ka(:), w(:)
  integer m, n, s, n1, n2, n3, n4, n5, n6
  double precision :: theta, usedummy
  usedummy = v(1) + va(1); usedummy = usedummy
  m = testcase % integers(1)
  n = testcase % integers(2)
  s = testcase % integers(3)
  n1 = n
  n2 = n1 + n*m
  n3 = n2 + n*n
  n4 = n3 + n*n*m
  n5 = n4 + s*n
  n6 = n5 + s
  allocate(w(n6), k(m), f0(n), F(n,m), A0(n,n), A(n,n,m), C(s,n), d(s), ka(m))
  w = testcase % reals
  k = u(1:m)
  f0 =         w(   1:n1)
  F  = reshape(w(n1+1:n2), [n,m])
  A0 = reshape(w(n2+1:n3), [n,n])
  A  = reshape(w(n3+1:n4), [n,n,m])
  C  = reshape(w(n4+1:n5), [s,n])
  d  =         w(n5+1:n6)
  call demo1_sumsq(m, n, s, k, f0, F, A0, A, C, d, theta, ka)
  ua(1:m) = ka
end subroutine F_rmd

PROGRAM TEST_DEMO1
  use rmd_dtesttools
  implicit none
  type(fdata) :: testcase
  double precision, allocatable :: u(:), k(:), f0(:), FF(:,:), A0(:,:), A(:,:,:),&
    &                              C(:,:) , d(:), w(:)
  double precision tol
  integer i, m, n, s, n1, n2, n3, n4, n5, n6
  external :: demo1_sumsq
  procedure(F_interf) :: F
  procedure(F_rmd_interf) :: F_rmd
  tol = 1.0d-6
  do m = 1,3
    do n = 2,5
      do s = 1,n-1
        allocate(u(m), k(m), f0(n), FF(n,m), A0(n,n), A(n,n,m), C(s,n), d(s))
        call random_number(k)
        call random_number(f0)
        call random_number(FF)
        call random_number(A0)
        call random_number(A)
        call random_number(C)
        call random_number(d)
        do i=1,n
          A0(i,i) = A0(i,i) + n
        enddo
        n1 = n
        n2 = n1 + n*m
        n3 = n2 + n*n;
        n4 = n3 + n*n*m
        n5 = n4 + s*n
        n6 = n5 + s
        allocate(w(n6))
        u = k
        w(   1:n1) = f0
        w(n1+1:n2) = reshape(FF, [n*m])
        w(n2+1:n3) = reshape(A0, [n*n])
        w(n3+1:n4) = reshape(A, [n*n*m])
        w(n4+1:n5) = reshape(C, [s*n])
        w(n5+1:n6) = d
        testcase = fdata(w, [m, n, s], ['d']) ! 'd' for dummy
        call rmd_dtestf(F, F_rmd, u, m, 1, tol, testcase)
        deallocate(w)
        deallocate(u, k, f0, FF, A0, A, C, d)
      enddo
    enddo
  enddo
END PROGRAM TEST_DEMO1
