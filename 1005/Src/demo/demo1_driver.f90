! DEMO1_DRIVER:  Driver program for the groundwater demo (example 1)
! 
! USAGE: 
!   make
!   ./demo1_driver
!
! Computes and displays the groundwater example sum of squares at two points, k1
! = (2,2) and k2 = (0.79, 1.83). These are the starting point and the
! approximate solution from the Matlab program demo_run_groundwater (in the
! BLAS-RMD matlab directory), which minimizes the sum of squares.

PROGRAM DEMO1_DRIVER
  ! use dispmodule
  use printutil
  implicit none
  integer, parameter :: m = 2, n = 6, s = 3
  double precision :: A0(n,n), A(n,n,m), h(7), f0(n), g(n), F(n,m), C(s,n),&
    &                 d(s), k1(m), k2(m), k1a(m), k2a(m), theta1, theta2
  character(10) fmt
  ! call disp_set(digmax=5, zeroas='0', advance='double')
  A0 = 0
  A(:,:,1) = reshape([ &
    &  4.,  -1.,   0., -2.,  0.,  0., &
    &  -1.,  2.,   0.,  0., -1.,  0., &
    &  0.,   0.,   0.,  0.,  0.,  0., &
    &  -1.,  0.,   0.,  4., -1.,  0., &
    &  0., -0.5,   0., -1.,  2.,  0., &
    &  0.,   0.,   0.,  0.,  0.,  0.], [6,6], order=[2,1])
  A(:,:,2) = reshape([ &
    &  0.,   0.,   0.,  0.,  0.,  0., &
    &  0.,   2.,  -1.,  0., -1.,  0., &
    &  0.,  -1.,   4.,  0.,  0., -2., &
    &  0.,   0.,   0.,  0.,  0.,  0., &
    &  0., -0.5,   0.,  0.,  2., -1., &
    &  0.,   0.,  -1.,  0., -1.,  4.], [6,6], order=[2,1])
  h = [100, 100, 75, 50, 25, 0, 0]
  f0 = 10
  g = [h(1), 0d0, h(7), h(2)+h(3), h(4), h(5)+h(6)]
  F = reshape([ &
    & g(1), g(2)/2, 0d0 , g(4), g(5)/2, 0d0 , &
    & 0d0 , g(2)/2, g(3), 0d0 , g(5)/2, g(6)], shape(F))
  C = reshape([ &
    & 0.375, 0.125, 0.0  , 0.375, 0.125, 0.0    , &
    & 0.0  , 0.250, 0.250, 0.0  , 0.250, 0.250  , &
    & 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.3333], shape(C), order=[2,1])
  d = [70d0, 35d0, 25d0 - h(5)/3 - h(6)/3]
  k1 = [2, 2]
  k2 = [0.7920, 1.8260]
  call demo1_sumsq(m, n, s, k1, f0, F, A0, A, C, d, theta1, k1a)
  call demo1_sumsq(m, n, s, k2, f0, F, A0, A, C, d, theta2, k2a)
  print '(A/)', 'BLAS-RMD example 1: Simplified groundwater model:'
  fmt = '(A,F0.4/)' 
  call printvec('k1 = ', k1)
  print fmt, 'theta1 = ', theta1
  call printvec('k1a = ', k1a)
  call printvec('k2 = ', k2)
  print fmt, 'theta2 = ', theta2
  call printvec('k2a = ', k2a)
  ! call disp(advance='yes')
  ! call disp('BLAS-RMD example 1: Simplified groundwater model:')
  ! call disp('k1 = '    , k1,     fmt = 'F0.3')
  ! call disp('theta1 = ', theta1, fmt = 'F0.3')
  ! call disp('k1a = '   , k1a,    fmt = 'F0.3')
  ! call disp('k2 = '    , k2,     fmt = 'F0.3')
  ! call disp('theta2 = ', theta2, fmt = 'F0.3')
  ! call disp('k2a = '   , k2a,    fmt = 'F0.3')
END PROGRAM DEMO1_DRIVER
