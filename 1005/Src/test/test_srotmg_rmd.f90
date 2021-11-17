PROGRAM TEST_SROTMG_RMD
  ! PURPOSE: Check that rotmg_rmd computes the correct adjoints.
  !
  ! USE:     test_srotmg_rmd
  !          test_srotmg_rmd -s
  !
  ! The -s flag is used to print out a table of the srotmg flag values and those
  ! of the scalings gamma that are tested.
  !
  ! NOTES: srotmg_rmd has several if- and while-statements, and the execution
  ! path depends on the values of flag, as well as the scaling done in srotmg.
  ! Thus the testing of srotmg_rmd includes several tests, to ensure full code
  ! coverage, and that several execution paths are tested.
  !
  ! It turns out that one must be careful in computing the numerical derivatives
  ! to compare with:
  !   1) It is necessary to choose the values where the adjoints are computed so
  !      that they are not too close to points where rotmg is not differentiable.
  !   2) An accurate finite difference stencil is used: a four-point stencil
  !      transpires to be sufficient
  !   3) The discretization step is chosen carefully
  !   4) Double precision is used to evaluate the numerical derivatives which
  !      are compared with the single precision analytic adjoint of srotmg_rmd.
  
  use rmd_stesttools
  implicit none
  procedure(F_interf)     :: F
  procedure(F_rmd_interf) :: F_rmd
  procedure(F_dp_interf)  :: F_dp
  integer :: t, i, dimu, dimv, flag
  type(fdata) :: testcase
  real :: mini, small, big, maxi, tol, tol0, u(4), v(7), step4p(4), step2p(4)
  real, allocatable :: testu(:,:)
  logical :: isdoubleversion

  isdoubleversion = .false.

  if (.not.isdoubleversion) then
    mini = 0
    small = 1e-7
    big = 1e7
    maxi = 0
    testu = reshape([ &
      2.,     2.,     2.,   3.,  &
      2.,     2.,     3.,   2.,  &
      small,  small,  2.,   3.,  &
      small,  small,  3.,   2.,  &
      small,  2.,     2.,   2.,  &
      small,  2.,     2e4,  2.,  &
      2.,     small,  2.,   2e4, &
      2.,     small,  2.,   2.,  &
      big,    big,    2.,   3.,  &
      big,    big,    3.,   2.], &
      [4, 10])
    tol0 = 1e-4
  else
    ! More comprehensive test
    mini = 1e-16
    small = 1e-9
    big = 1e9
    maxi = 1e17
    testu = reshape([ &
      2.,     2.,     2.,   3.,  &
      2.,     2.,     3.,   2.,  &
      small,  small,  2.,   3.,  &
      small,  small,  3.,   2.,  &
      small,  2.,     2.,   2.,  &
      small,  2.,     2e4,  2.,  &
      small,  mini,   3.,   2.,  &
      mini,   mini,   1.,   3.,  &
      mini,   mini,   3.,   1.,  &
      2.,     small,  2.,   2e4, &
      2.,     small,  2.,   2.,  &
      big,    big,    2.,   3.,  &
      big,    big,    3.,   2.,  &
      maxi,   big,    2.,   3.], &
      [4, 14])
    tol0 = 1e-5
  endif

  if (sflagoncommandline()) print "('Flag     Gamma1         Gamma2')"
  dimu = 4
  do t = 1, size(testu, 2)
    u = testu(:, t)
    step2p = 5e-3 ! 2-point stencil
    where(u == small) step2p = small/500
    where(u == big)   step2p = big/500
    step4p = 5e-2 ! 4-point stencil
    where(u == mini)  step4p = mini/50
    where(u == small) step4p = small/50
    where(u == big)   step4p = big/50
    where(u == maxi)  step4p = maxi/50
    tol = tol0
    if (sflagoncommandline()) call scaling_table(u)
    ! Use testcase % integers to transfer flag from rotmg to rotmg_rmd
    testcase = fdata([real::], [0], [character::])
    do i = 1,3
      call F(testcase, u, v)
      flag = testcase % integers(1)
      if (flag <  0) dimv = 7
      if (flag >= 0) dimv = 5
      if (isdoubleversion) then
        call rmd_stestf(F, F_rmd, u,  dimu, dimv, tol, testcase, step = step4p, &
          &             fourpoint = [.true., .true., .true., .true.])
      else
        call rmd_stestf(F, F_rmd, u,  dimu, dimv, tol, testcase, step = step2p, &
          &             F_dp = F_dp)
      endif
      call random_number(u)
      u = testu(:,t)*(1 + (u - 0.5)/4)
      tol = tol0*20
    enddo
  enddo

contains
  function sflagoncommandline() result(sf)
    character(3) cmdarg
    logical :: sf
    sf = .false.
    if (command_argument_count() > 0) then
      call get_command_argument(1, cmdarg)
      if (cmdarg == "-s") then
        sf = .true.
      else
        stop "use -s to display table of tested flags and scalings"
      endif
    endif
  end function sflagoncommandline

  subroutine scaling_table(u)
    ! Print table of combinations of Flag and scalings that are tested
    ! (used to ascertain that all code segments in rotmg_rmd are actually tested)
    integer :: flag
    real :: u(*), d(2), x1, y1, param(5), da(2), x1a, y1a, parama(8), gam1, gam2
    d = u(1:2)
    x1 = u(3)
    y1 = u(4)
    call srotmg(d(1), d(2), x1, y1, param)
    da = [1.0, 1.0]
    x1a = 1.0
    parama(1:5) = [2.0, 1.0, 1.0, 1.0, 1.0]
    call srotmg_rmd(d(1), d(2), x1, param, da(1), da(2), x1a, y1a, parama)
    flag = nint(parama(6))
    gam1 = parama(7)
    gam2 = parama(8)
    print '(I0,T7,1PG14.7,1X,G14.7)', flag, gam1, gam2
  end subroutine scaling_table
END PROGRAM TEST_SROTMG_RMD

SUBROUTINE F(testcase, u, v)
  use rmd_stesttools
  implicit none
  procedure(set_KJ_interf) :: set_KJ
  type(fdata), intent(inout) :: testcase
  real, intent(in)  :: u(*)
  real, intent(out) :: v(*)
  real :: d1, d2, x1, y1, param(5)
  integer :: flag
  integer, allocatable :: K(:), J(:)
  !
  d1 = u(1)
  d2 = u(2)
  x1 = u(3)
  y1 = u(4)
  call srotmg(d1, d2, x1, y1, param)
  v(1) = d1
  v(2) = d2
  v(3) = x1
  flag = nint(param(1))
  call set_KJ(flag, K, J)
  v(J) = param(K)
  testcase % integers(1) = flag
END SUBROUTINE F

SUBROUTINE F_DP(testcase, u, v)
  use rmd_stesttools
  implicit none
  procedure(set_KJ_interf) :: set_KJ
  type(fdata), intent(inout) :: testcase
  double precision, intent(in)  :: u(*)
  double precision, intent(out) :: v(*)
  double precision :: d1, d2, x1, y1, param(5)
  integer :: flag
  integer, allocatable :: K(:), J(:)
  !
  d1 = u(1)
  d2 = u(2)
  x1 = u(3)
  y1 = u(4)
  call drotmg(d1, d2, x1, y1, param)
  v(1) = d1
  v(2) = d2
  v(3) = x1
  flag = nint(param(1))
  call set_KJ(flag, K, J)
  v(J) = param(K)
  testcase % integers(1) = flag
END SUBROUTINE F_DP

SUBROUTINE F_RMD(testcase, u, v, ua, va)
  use rmd_stesttools
  implicit none
  procedure(set_KJ_interf) :: set_KJ
  type(fdata), intent(inout) :: testcase
  real, intent(in)  :: u(*), v(*), va(*)
  real, intent(inout) :: ua(*)
  real :: d1, d2, x1, param(5), parama(5)
  real :: d1a, d2a, x1a, y1a
  integer :: flag
  integer, allocatable :: K(:), J(:)
  !
  flag = testcase % integers(1)
  d1 = v(1)
  d2 = v(2)
  x1 = v(3)
  d1a = va(1)
  d2a = va(2)
  x1a = va(3)
  call set_KJ(flag, K, J)
  param(1) = flag
  param(K) = v(J)
  parama(1) = flag
  parama(K) = va(J)
  call srotmg_rmd(d1, d2, x1, param, d1a, d2a, x1a, y1a, parama)
  ua(1) = d1a
  ua(2) = d2a
  ua(3) = x1a
  ua(4) = y1a
END SUBROUTINE F_RMD

subroutine set_KJ(flag, K, J)
  integer, intent(in) :: flag
  integer, intent(out), allocatable :: K(:), J(:)
  select case(flag)
  case(0)
    K = [3,4]
    J = [4,5]
  case(1)
    K = [2,5]
    J = [4,5]
  case(-1)
    K = [2, 3, 4, 5]
    J = [4, 5, 6, 7]
  case(-2)
    K = [integer::]
    J = [integer::] 
  case default
    stop 'Unexpected flag'
  end select
end subroutine set_KJ
