! DLYAP  Solve discrete Lyapunov equation via SLICOT
! 
! CALL DLYAP(transp, A, Sig, S) solves the discrete Lyapunov equation
! 
!             S - op(A)*S*op(A)' = Sig,
!             
! for S, where op(A) is A or A' depending on whether transp is 'n' or 't'. It
! calls the SLICOT subroutine SB03MD.

subroutine dlyap(transp, n, A, Sig, S)
  implicit none
  integer, parameter :: dble = kind(0d0)
  character, intent(in)     :: transp   ! work with A or A'
  integer, intent(in)       :: n        ! matrix size
  real(dble), intent(in)    :: A(n,n)   ! n by n, coefficients
  real(dble), intent(inout) :: Sig(n,n) ! n by n, rhs, symmetric, stored in lower triangle
  real(dble), intent(inout) :: S(n,n)   ! n by n, solution, symmetric, both triangles stored

  character(1) :: dico, job, fact, trana
  integer lda, ldu, ldc, ldwork, info
  real(dble) :: scale, sep, ferr
  real(dble), allocatable :: U(:,:), wr(:), wi(:), dwork(:), AA(:,:)
  integer, allocatable :: iwork(:)
  
  dico = 'D'
  job = 'X'
  fact = 'N'
  if (transp == 'n' .or. transp == 'N') then ! opposite to dlyap
    trana = 'T'
  else
    trana = 'N'
  endif
  lda = n
  ldu = n
  ldc = n
  ldwork = 2*n*(n+2)
  allocate (U(n,n), wr(n), wi(n), iwork(n*n), dwork(ldwork), AA(n,n))
  call dcopy(n*n, Sig, 1, S, 1)
  call dcopy(n*n, A, 1, AA, 1)
  call sb03md(dico, job, fact, trana, n, S, ldc, AA, lda, U, ldu, scale &
    ,         sep, ferr, wr, wi, iwork, dwork, ldwork, info)
  if (info < 0) then
    print '(A,I0,A)', 'On entry to sb03md argument number ', -info &
      ,              'had an illegal value'
    stop
  elseif (info > 0) then
    print '(A,I0)', 'sb03md failed to solve the discrete Lyapunov '// &
      'equation, exit with info = ', info
    stop
  endif
  call dscal(n*n, -1/scale, S, 1)
end subroutine dlyap
