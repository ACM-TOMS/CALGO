subroutine demo1_sumsq(m, n, s, k, f0, F, A0, A, C, d, theta, ka)
  !use dispmodule
  implicit none
  character(2) :: select
  integer, intent(in) :: m, n, s
  double precision, intent(in) :: k(m), f0(n), F(n,m), A0(n,n), A(n,n,m),&
    &                             C(s,n), d(s)
  double precision, intent(out) :: theta, ka(m)
  integer :: j, info
  double precision :: thetaa, dummy(1)
  double precision, external :: ddot
  integer, allocatable :: ipiv(:)
  double precision, allocatable :: g(:),  r(:),  u(:),  B(:,:), z(:),&
    &                              ga(:), ra(:), ua(:), Ba(:,:)
  allocate(g(m), ga(m), r(s), ra(s), u(n), ua(n), B(n,n), Ba(n,n), z(n), ipiv(n))
  ! COMPUTE SUM-OF-SQUARES
  B = A0;
  do j=1,m
    call daxpy(n*n, k(j), A(1,1,j), 1, B, 1)         ! B = B + k(j)*A(:,:,j);
  enddo
  g = f0;
  call dgemv('n', n, m, 1d0, F, n, k, 1, 1d0, g, 1)  ! g := F*k + g
  call dgetrf(n, n, B, n, ipiv, info)            ! linsys: PLU-factorize B
  if (info /= 0) stop 'singular matrix'          ! linsys
  u = g                                          ! linsys
  call dgetrs('n', n, 1, B, n, ipiv, u, n, info) ! linsys: Solve B*u = g
  r = d;
  call dgemv('n', s, n, 1d0, C, s, u, 1, -1d0, r, 1) ! r := C*u - d
  theta = ddot(s, r, 1, r, 1);
  
  ! COMPUTE ADJOINT OF k
  thetaa = 1; ra = 0; ua = 0; Ba = 0; ka = 0
  select = '11'
  call ddot_rmd(s, r, 1, r, 1, thetaa, ra, ra, select)
  call dgemv_rmd('n', s, n, 1d0, C, s, u, 1, -1d0, 1, dummy, ua, ra, '010')
  z = ua
  call dgetrs('t', n, 1, B, n, ipiv, z, n, info)  ! linsys_rmd: Solve B'*z = ua
  ga = z                                          ! linsys_rmd
  call dger(n, n, -1d0, z, 1, u, 1, Ba, n)        ! linsys_rmd: Ba := -z*u'
  call dgemv_rmd('n', n, m, 1d0, F, n, k, 1, 1d0, 1, dummy, ka, ga, '010')
  do j=m,1,-1
    ka(j) = ka(j) + ddot(n*n, Ba, 1, A(1,1,j), 1)
  enddo
end subroutine demo1_sumsq
