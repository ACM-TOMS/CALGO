function demo2_likeli (n, r, Sig, A, X, Siga, Aa, compute_adjoint) result(l)
  !use dispmodule
  use printutil
  implicit none
  integer, intent(in)             :: n, r
  double precision, intent(in)    :: A(r,r), X(r,n)
  double precision, intent(inout) :: Sig(r,r)
  double precision, intent(out)   :: Siga(r,r), Aa(r,r)
  logical, intent(in)             :: compute_adjoint
  integer                         :: infochol, i
  double precision                :: v, sumdiag, l, ddot, la, va, sumdiaga,&
    &                                wrk(1) =[0.0d0], dummy(1)
  double precision, parameter     :: pi = 4d0*atan(1d0)
  double precision, allocatable   :: S(:,:), LS(:,:), LSig(:,:), u(:), ua(:), &
    &                                W(:,:), Wa(:,:), Sa(:,:)
  allocate(S(r,r), LS(r,r), LSig(r,r), u(r), ua(r), W(r,n-1), Wa(r,n), Sa(r,r))
  ! Total memory 4r^2 + 2rn

  ! CALCULATE LIKELIHOOD
  call sym(r, Sig, Sig)
  call dlyap('n', r, A, Sig, S)
  ! implicit: S := isym(S)
  LS = S
  call dpotrf('L', r, LS, r, infochol)
  LSig = Sig
  call dpotrf('L',r, LSig, r, infochol)
  u = X(:,1)
  W = X(:,2:n)
  call dgemm('n', 'n', r, n-1, r, -1d0, A, r, X, r, 1d0, W, r)
  call dtrsv('L', 'n', 'n', r, LS, r, u, 1)
  call dtrsm('L', 'L', 'n', 'n', r, n-1, 1d0, LSig, r, W, r)
  v = ddot(r, u, 1, u, 1) + ddot(r*(n-1), W, 1, W, 1)
  sumdiag = 0
  do i = 1,r
    sumdiag = sumdiag + log(LS(i,i)) + (n-1)*log(LSig(i,i))
  enddo
  l = -(n*r*log(2*pi) + 2*sumdiag + v)/2

  if (compute_adjoint) then
    la = 1d0
    va = -la/2
    sumdiaga = -la
    ! Use same variable for Siga and LSiga as well as for Sa and LSa
    Siga = 0; Sa = 0; Wa = 0; ua = 0; Aa = 0
    ! Note: if z = sum(f(x(j))) then xa(j) = za*f'(x(j))
    do i = r,1,-1
      Siga(i,i) = (n-1)*sumdiaga/LSig(i,i)
      Sa(i,i) = sumdiaga/LS(i,i)
    enddo
    call ddot_rmd(r*(n-1), W, 1, W, 1, va, Wa, Wa, '11')
    call ddot_rmd(r, u, 1, u, 1, va, ua, ua, '11')
    call dtrsm_rmd('L', 'L', 'n', 'n', r, n-1, 1d0, LSig, r, W, r, Siga, Wa, wrk, '11')
    call dtrsv_rmd('L', 'n', 'n', r, LS, r, u, 1, Sa, ua, wrk, '11')
    !sel = '100': Need Aa but not adjoint of X and W
    call dgemm_rmd('n', 'n', r, n-1, r, -1d0, A, r, X, r, 1d0, r, Aa, dummy, Wa, '100')
    call dpotrf_rmd('L', r, LS, r, Sa)
    call dpotrf_rmd('L', r, LSig, r, Siga)
    call isym_rmd(r, Sa)
    ! if S solves S - A*S*A' = Sig then Siga solves Siga - A'*Siga*A = Sa
    call dlyap('t', r, A, Sa, Sa) ! Overwrite Sa with Siga contribution from dlyap
    call dsymm('L', 'L', r, r, 1d0, Sa, r, A, r, 0d0, W, r) ! Use W to store Sa*A
    call dsymm('R', 'L', r, r, 2d0, S, r, W, r, 1d0, Aa, r) ! Aa += 2*Sa*A*S
    call sym_rmd(r, Sa)
    Siga = Siga + Sa
  endif
end function demo2_likeli

subroutine sym(n, L, S)
  ! Symmetric matrix from lower part, i.e. S:=sym(L). It is allowed to have L=S.
  integer, intent(in) :: n
  double precision, intent(inout) :: L(n,*)
  double precision, intent(inout) :: S(n,*)
  integer :: i
  do i = 1,n
    S(i,1:i) = L(i,1:i)
    S(1:i-1,i) = L(i,1:i-1)
  enddo
end subroutine sym

subroutine sym_rmd(n, S)
  ! Overwrites S with its adjoint for the sym operation, 
  integer, intent(in) :: n
  double precision, intent(inout) :: S(n,*)
  integer :: i
  do i=1,n
    S(i,i+1:n) = 0d0
    call dscal(n - i, 2d0, S(i+1:n,i), 1)
  enddo
end subroutine sym_rmd

subroutine isym_rmd(n, S)
  ! Overwrites S with its adjoint for the isym operation, 
  integer, intent(in) :: n
  double precision, intent(inout) :: S(n,*)
  integer :: i
  do i=1,n
    call dscal(n-i, 0.5d0, S(i+1:n, i), 1)
    S(i, i+1:n) = S(i+1:n, i)
  enddo
end subroutine isym_rmd
