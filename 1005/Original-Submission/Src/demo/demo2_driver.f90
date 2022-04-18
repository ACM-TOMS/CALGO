! DEMO2_DRIVER:  Driver program for VAR-1 likelihood demo (example 2).
! 
! The program computes the likelihood function and optionally its adjoint for a
! given set of parameters (the default is to compute the adjoint). Any of the
! following command syntaxes may be used to run it:
!
!      $ ./demo2_driver <datafile>
!      $ ./demo2_driver <r> <n>
!      $ ./demo2_driver -n <datafile>
!      $ ./demo2_driver -n <r> <n>
!
! where
!      <datafile> is generated by the program demo_var_init.m or contains data
!      of a vector autoregressive p=1 (VAR-1) time series, in a format as is
!      generated by demo_var_init.m
!
!      <r> is the dimension of the vector at each timestep and <n> is the
!      number of terms in the series. With this syntax data for A is
!      selected randomly, X(i,j) is set to (i+j)/(n+m) and Sigma is set to
!      a Lehmer matrix
!
!      The switch -n may be used to suppress the adjoint computation.
!
! The demo directory contains four example datafiles, named var_<r>_<n>.dat,
! where r is the dimension of the series and n its length.

PROGRAM DEMO2_DRIVER
  !use dispmodule
  use printutil
  implicit none
  integer n, r, i, j, cac, offset
  logical compute_adjoint
  double precision, allocatable :: A(:,:), Sigma(:,:), X(:,:), Aa(:,:), Sigmaa(:,:)
  double precision lval
  double precision, external :: demo2_likeli
  character(30) :: filename, rs, ns, arg1
  !
  !call disp_set(zeroas = '0', digmax = 5, advance = 'double')
  cac = command_argument_count()
  offset = 0
  if (cac == 0) then
    print *,'Usage examples:'
    print *,'  demo2_driver var_2_3.dat     # Read parameters from file'
    print *,'  demo2_driver 10 1000         # Generate parameters for r=10, n=1000'
    print *,'  demo2_driver -n var_2_3.dat  # No adjoint, parameters from file'
    print *,'  demo2_driver -n 10 1000      # No adjoint, generate parameters'
    stop
  end if

  call get_command_argument(1, arg1)
  if (len(arg1) == 0) stop 'zero-length argument'
  if (arg1 == '-n') then
    compute_adjoint = .false.
    offset = offset + 1
  else
    compute_adjoint = .true.
  endif
  if (cac - offset == 1) then
    call get_command_argument(offset + 1, filename)
    open(unit=10, file=filename, status='old')
    read(10,*) r
    read(10,*) n
  else
    call get_command_argument(offset + 1, rs)
    call get_command_argument(offset + 2, ns)
    read(rs,*) r
    read(ns,*) n
  endif
  if (n < r) stop 'Cannot have n < r'
  allocate(A(r,r), Sigma(r,r), X(r,n), Aa(r,r), Sigmaa(r,r))
  if (cac - offset == 1) then
    read(10,*) ((X(i,j), j=1,n), i=1,r)
    read(10,*) ((Sigma(i,j), j=1,r), i=1,r)
    read(10,*) ((A(i,j), j=1,r), i=1,r)
  else
    call setx(X)
    call lehmer(Sigma)
    call random_number(A)
    call dscal(r*r, 0.1d0/r, A, 1) 
  endif
  !call disp(advance='yes')
  ! call disp('BLAS-RMD example 2: Vector-autoregressive (VAR) time series likelihood:')
  ! call disp('A = ', A)
  ! call disp('X = ', X)
  ! call disp('Sigma = ', Sigma)
  print '(A/)','BLAS-RMD example 2: Vector-autoregressive (VAR) time series like&
    &lihood:'
  if (r <= 5) then
    call printmat('A = ', A)
    call printmat('X^T = ', transpose(X))
    call printmat('Sigma = ', Sigma)
  else
    print '(A)', 'First 5 by 5 parts of input matrices:'
    call printmat('A = ', A(1:5, 1:5))
    call printmat('X = ', X(1:5, 1:5))
    call printmat('Sigma = ', Sigma(1:5, 1:5))
  end if
  !
  lval = demo2_likeli(n, r, Sigma, A, X, Sigmaa, Aa, compute_adjoint)
  !
  print '(A, F0.4/)', 'likelihood = ', lval
  ! call disp('likelihood = ', lval)
  ! call disp('Sigmaa = ', Sigmaa)
  ! call disp('Aa = ', Aa)
  if (compute_adjoint) then
    if (r <= 5) then
      call printmat('Sigma-adjoint = ', Sigmaa)
      call printmat('A-adjoint = ', Aa)
    else
      print '(A)', 'First 5 by 5 parts of adjoint matrices:'
      call printmat('Sigmaa = ', Sigmaa(1:5,1:5))
      call printmat('Aa = ', Aa(1:5,1:5))
    endif
  endif
  
CONTAINS
  subroutine lehmer(A)
    double precision, intent(out) :: A(:,:)
    integer :: i, j, n
    n = size(A,1)
    do i = 1,n
      do j = 1,i
        A(i,j) = dble(j)/i
        A(j,i) = A(i,j)
      enddo
    enddo
  end subroutine lehmer

  subroutine setx(X)
    double precision, intent(out) :: X(:,:)
    integer :: i, j, m, n
    m = size(X,1)
    n = size(X,2)
    do i = 1,m
      do j = 1,n
        X(i,j) = dble(i + j)/(m + n)
      enddo
    enddo
  end subroutine setx
END PROGRAM DEMO2_DRIVER