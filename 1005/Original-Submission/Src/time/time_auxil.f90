SUBROUTINE RUNBLAS(blasrout, A, X, Y, n)
  ! Run a blas routine
  ! PARAMETERS
  use param
  use tictoc
  character(5), intent(in)    :: blasrout ! e.g. dgemv
  real(dble),   intent(in)    :: A(*)     ! lhs BLAS-matrix operand
  real(dble),   intent(in)    :: X(*)     ! rhs BLAS-operand (vector or matrix)
  real(dble),   intent(inout) :: Y(*)     ! BLAS-result (vector or matrix)
  integer,      intent(in)    :: n        ! BLAS matrices are n by n
  real(dble) :: alpha, beta
  integer :: bandw
  !
  alpha = 1.0000001d0
  beta  = 0.9999999d0 ! so close to 1 to prevent underflow/overflow
  select case(blasrout)
  case('dgemv')
    call dgemv('n', n, n, alpha, A, n, X, 1, beta, Y, 1)
  case('dgemm')
    call dgemm('n', 'n', n, n, n, alpha, A, n, X, n, beta, Y, n)
  case('dsyrk')
    call dsyrk('u', 'n', n, n, alpha, A, n, beta, X, n)
  case('dtrsm')
    Y(1:n**2) = X(1:n**2) ! necessary to prevent underflow/overflow
    call dtrsm('l', 'u', 'n', 'u', n, n, alpha, A, n, Y, n)    
  case('dtrmv')
    call dtrmv('l', 'n', 'n', n, A, n, X, 1)
  case('dgbmv')
    bandw = n/3
    call dgbmv('n', n, n, bandw, bandw, alpha, A, n, X, 1, beta, Y, 1)
  case default
    stop 'Unknown blas function'
  end select
END SUBROUTINE RUNBLAS

SUBROUTINE RUNBLASRMD(blasrout, A, X, Y, Aa, Xa, Ya, n)
  ! Run a BLAS-RMD routine
  ! PARAMETERS
  use param
  character(5), intent(in)    :: blasrout ! e.g. dgemv
  real(dble),   intent(in)    :: A(*)     ! lhs BLAS-matrix operand
  real(dble),   intent(in)    :: X(*)     ! rhs BLAS-operand (vector or matrix)
  real(dble),   intent(inout) :: Y(*)     ! BLAS-result (vector or matrix)
  real(dble),   intent(inout) :: Aa(*)    ! A-adjoint
  real(dble),   intent(inout) :: Xa(*)    ! X-adjoint
  real(dble),   intent(inout) :: Ya(*)    ! Y-adjoint
  integer,      intent(in)    :: n        ! BLAS matrices are n by n
  real(dble) :: alpha, beta
  integer :: bandw
  !
  real(dble) useY
  useY = Y(1); useY = useY
  !
  alpha = 1.0000001d0
  beta  = 0.9999999d0 ! so close to 1 to prevent underflow/overflow
  select case(blasrout)
  case('dgemv')
    call dgemv_rmd('n', n, n, alpha, A, n, X, 1, beta, 1, Aa, Xa, Ya, '111')
  case('dgemm')
    call dgemm_rmd('n', 'n', n, n, n, alpha, A, n, X, n, beta, n, Aa, Xa, Ya, '111')
  case('dsyrk')
    call dsyrk_rmd('u', 'n', n, n, alpha, A, n, beta, n, Aa, Xa, '111')
  case('dtrsm')
    ! use Y as workspace
    call dtrsm_rmd('l', 'u', 'n', 'u', n, n, alpha, A, n, X, n, Aa, Ya, Y, '11')
  case('dtrmv')
    call dtrmv_rmd('l', 'n', 'n', n, A, n, X, 1, Aa, xa, '11')
  case('dgbmv')
    bandw = n/3
    call dgbmv_rmd('n', n, n, bandw, bandw, alpha, A, n, X,1,beta,1, Aa,Xa,Ya,'111')
  case default
    stop 'Unknown blas function'
  end select
END SUBROUTINE RUNBLASRMD
