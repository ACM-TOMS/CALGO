MODULE PARAM
  integer, parameter :: &
    &  dble = kind(0d0), &
    &  byt4 = selected_int_kind(9), &
    &  byt8 = selected_int_kind(15), &
    &  bytx = (byt4 + 1 + 2*byt8 + sign(byt4 + 1, -byt8))/2
  ! Kinds are positive or -1 if not possible, so bytx will signify 8-byte integers if
  ! these exist, otherwise 4-byte integers.
END MODULE PARAM

MODULE TICTOC
  use param
  private
  public tic, toc
  integer(bytx) :: count, rate, count1
  
CONTAINS
  ! function sep000(n) result(s)
  !   ! Tostring with . as thousands separator
  !   character(len=:), allocatable :: s
  !   integer(byt8), intent(in) :: n
  !   integer(byt8) :: nv(8), nn, i
  !   nn = n
  !   i = 8
  !   do while(nn > 0)
  !     nv(i) = mod(nn,1000)
  !     i = i-1
  !     nn = nn/1000
  !   enddo
  !   call tostring_set(sep = '.')
  !   s = tostring(nv(i+1)) // '.' // tostring(nv(i+2:8), 'I3.3')
  ! end function sep000
  
  subroutine tic()
    implicit none
    call system_clock(count, rate)
  end subroutine tic

  function toc() result(toctime)
    implicit none
    real(dble) toctime
    call system_clock(count1, rate)
    toctime = real(count1 - count, dble)/rate
  end function toc
END MODULE TICTOC

MODULE TIMING
  use tictoc
  use param
  private
  public initw, timing_info, timeblas

  type timing_info
    ! each element has [time in seconds, number of calls]
    real(dble)          &
      &   one_blas(2),  &
      &   many_blas(2), &
      &   one_rmd(2),   &
      &   many_rmd(2)
  end type timing_info

CONTAINS
  subroutine initw(w, nw)
    integer, intent(in) :: nw
    real(dble), intent(inout) :: w(nw)
    call random_number(w)
    !w = 0.1
  end subroutine initw

  subroutine run(rmd, blasrout, w, nw, iw, n, skip)
    ! Run a blas or BLAS-RMD routine
    ! PARAMETERS:
    logical, intent(in)       :: rmd      ! run BLAS-RMD if true, otherwise blas
    character(5), intent(in)  :: blasrout ! e.g. 'dgemv' to run dgemv/dgemv_rmd
    integer, intent(in)       :: n        ! BLAS matrices are n by n
    integer, intent(in)       :: iw       ! pointer into current data group in w
    integer, intent(in)       :: skip(:)  ! number of elements in A, X and Y
    integer, intent(in)       :: nw       ! size of w
    real(dble), intent(inout) :: w(nw)    ! workspace with blas matrix/vector arguments
    !
    integer :: iA, iX, iY, iAa, iXa, iYa
    iA  = iw
    iX  = iA + skip(1)
    iY  = iX + skip(2)
    if (rmd) then ! time BLAS-RMD
      iAa = iY + skip(3)
      iXa = iAa + skip(1)
      iYa = iXa + skip(2)
      call runblasrmd(blasrout, w(iA), w(iX), w(iY), w(iAa), w(iXa), w(iYa), n)
    else
      call runblas(blasrout, w(iA), w(iX), w(iY), n)
    endif
  end subroutine run

  function time1(type, blasrout, w, nw, n, skip, wskip, allottedtime) result(time)
    ! Time calls of a blas or BLAS-RMD routine repeatedly (for improved accuracy).
    ! Returns time for each call in seconds
    ! PARAMETERS:
    character(*), intent(in)    :: type         ! 'blas' or 'rmd'
    character(*), intent(in)    :: blasrout     ! 'dgemv'/'dgemm'/...
    integer,      intent(in)    :: nw           ! size of w
    real(dble),   intent(inout) :: w(nw)        ! Work space with BLAS-matrices/vectors
    integer,      intent(in)    :: n            ! BLAS matrices are n by n
    integer,      intent(in)    :: skip(:)      ! stride between matrices/vectors
    integer,      intent(in)    :: wskip        ! stride between groups of data cases in w
    real(dble),   intent(in)    :: allottedtime ! Run repeatedly for at least this time
    real(dble)                  :: time(2)      ! runtimes for this case
    !
    integer :: r, reps, atleast, additional, iw
    logical :: rmd, firsttime
    real(dble) :: toct
    !
    select case(type)
    case('blas')
       rmd = .false.
    case('rmd')
       rmd = .true.
    end select
    r = 0
    reps = 0
    additional = 1
    iw = 1
    call tic()
    firsttime = .true.
    do while(.true.)
      do while (r < reps + additional)
        call run(rmd, blasrout, w, nw, iw, n, skip)
        r = r + 1
        if (firsttime) then
          firsttime = .false.
          toct = toc()
          if (toct > 0.2 .and. toct > 0.1*allottedtime) then
            r = 0
            call tic()
            cycle
          endif
        endif
        iw = iw + wskip
        if (iw + wskip - 1 > nw) iw = 1
      end do
      reps = reps + additional
      time(1) = toc()
      atleast = max(reps + 1, nint(reps*1.1))
      additional = nint(max(reps*allottedtime/time(1), real(atleast,dble))) - reps;
      if (time(1) > allottedtime*0.9) exit
      if (time(1) < 1d-6) then
        additional = 100*reps;
      elseif (time(1) < 20d-6) then
        additional = max(additional, 5*reps)
      endif
    end do
    time(1) = time(1)/r
    time(2) = r
  end function time1

  function timeblas(rout, w, nw, n, skip, wskip, allottedtime) result(time)
    ! Time four or six blasroutines (dgemv, dgemm, dsyrk, dtrsm, and optionally,
    ! dtrmv, dgbmv) along with their rmd counterparts
    character(*), intent(in)    :: rout         ! dgemv/dgemm/dsyrk/dtrsm/dtrmv/dgbmv
    integer,      intent(in)    :: nw           ! size of w
    real(dble),   intent(inout) :: w(nw)        ! Work space with BLAS-matrices/vectors
    integer,      intent(in)    :: n            ! BLAS matrices are n by n
    integer,      intent(in)    :: skip(:)      ! stride between matrices/vectors
    integer,      intent(in)    :: wskip        ! stride between groups of data cases in w
    real(dble),   intent(in)    :: allottedtime ! Run repeatedly for at least this time
    type(timing_info)           :: time         ! runtimes for all six cases
    !
    time % one_blas    = time1('blas', rout, w, nw, n, skip, 0,     allottedtime)
    time % many_blas   = time1('blas', rout, w, nw, n, skip, wskip, allottedtime)
    time % one_rmd     = time1('rmd',  rout, w, nw, n, skip, 0,     allottedtime)
    time % many_rmd    = time1('rmd',  rout, w, nw, n, skip, wskip, allottedtime)
  end function timeblas
END MODULE TIMING

MODULE PRINTING
  use param
  use timing
CONTAINS
  function tostr(n) result(pstring)
    integer, intent(in) :: n(:)
    character(50) :: pstring, fmt
    fmt = '(   (I0,","),I0)'
    write(fmt(2:4),'(I3)') ubound(n,1)-1
    write(pstring, fmt) n
  end function tostr

  subroutine time2vec(time, t, nc)
    type(timing_info), intent(in) :: time
    real(dble), intent(out) :: t(4)
    integer, intent(out) :: nc(4)
    t(1) = time % one_blas(1)
    t(2) = time % many_blas(1)
    t(3) = time % one_rmd(1)
    t(4) = time % many_rmd(1)
    nc(1) = nint(time % one_blas(2))
    nc(2) = nint(time % many_blas(2))
    nc(3) = nint(time % one_rmd(2))
    nc(4) = nint(time % many_rmd(2))
  end subroutine time2vec
  
  function timesum(time) result(tsum)
    type(timing_info), intent(in) :: time
    real(dble) :: tsum
    real(dble) :: t(4)
    integer :: nc(4)
    call time2vec(time, t, nc)
    tsum = sum(t*nc)
  end function timesum

  subroutine print(col1, n, time, factor, nthr, nspread)
    character(*),      intent(in)           :: col1
    integer,           intent(in), optional :: n
    real(dble),        intent(in), optional :: factor
    type(timing_info), intent(in), optional :: time
    integer,           intent(in), optional :: nspread
    integer,           intent(in), optional :: nthr
    !
    real(dble) :: t(4)
    integer :: nc(4), ps_per_flop(4)
    character(*), parameter :: &
      & fmt(3) = &
      &   ['(   T11, 8X, 5X, 7X, A14, 2X, A14,  3X, A12,   3X, A)', &
      &    '(A, T11, A8, A5, A7, A14, 2X, A14,  3X, A12,   3X, A)', &
      &    '(A, T11, I8, I5, I7, 2I7, 2X, 2I7,  3X, 2F6.1, 3X, A)'], &
      !
      & h1(3) = ['       BLAS-  ', &
      &          '    subroutine', &
      &          '   same spread'], &
      !
      & h2(3) = [' Corresponding', &
      &          'RMD-subroutine', &
      &          '   same spread'], &
      !
      & h3(3) = ['  RMD-time /', &
      &          '   BLAS-time', &
      &          ' same spread'], &
      & h4(3) = ['               ', &
      &          '               ', &
      &          'number of calls']

    if(.not.present(n)) then
      print fmt(1), h1(1), h2(1), h3(1), h4(1)
      print fmt(1), h1(2), h2(2), h3(2), h4(2)
      print fmt(2), 'Operation', 'nspread', 'thr', 'n', h1(3), h2(3), h3(3), h4(3)
    else
      call time2vec(time, t, nc)
      ps_per_flop(1:2) = nint(t(1:2)*factor)
      ps_per_flop(3:4) = nint(t(3:4)*factor/2)
      print fmt(3), col1, nspread, nthr, n, ps_per_flop, t(3:4)/t(1:2), trim(tostr(nc))
    end if
  end subroutine print

  subroutine printhelp
    print *, 'TIMERMD'
    print *, '   Time BLAS routines and corresponding reverse mode derivatives (RMD)'
    print *
    print *, 'USAGE:'
    print *, '   timermd               print this help'
    print *, '   timermd -             use defaults for all parameters'
    print *, '   timermd -nN1,N2...    matrix sizes to time, default N1=100, N2=1000'
    print *, '   timermd -tNt1,Nt2...  number of threads, default 1'
    print *, '   timermd -sS           allot S sec for each timing, default 1'
    print *, '                         (the default, -n100,1000 counts as two timings)'
    print *, '   timermd -mMemsiz      use Memsiz elements for spread matrices, default'
    print *, '                         as needed for one spread for the maximum N'
    print *
    print *, 'EXAMPLES:'
    print *, '   timermd -s2 -m1.99G -n1000'
    print *, '   timermd -m800M -n100,1000,10000 -t1,2'
  end subroutine printhelp

  subroutine printsummary(memsiz, totaltime)
    integer, intent(in) :: memsiz
    real(dble), intent(in) :: totaltime
    real(dble) :: bytes
    bytes = memsiz*8d0
    if (bytes < 0) then
    elseif (bytes < 1024) then
      print '("Memory: ",I0," bytes. Total time: ",F0.1," sec.")', nint(bytes), totaltime
    elseif (bytes < 2**20) then
      print '("Memory: ",F0.1," KB. Total time: ",F0.1, " sec.")', bytes/2**10, totaltime
    elseif (bytes < 2**30) then
      print '("Memory: ",F0.1," MB. Total time: ",F0.1," sec.")', bytes/2**20, totaltime
    else
      print '("Memory: ",F0.1," GB. Total time: ",F0.1," sec.")', bytes/2**30, totaltime
    endif
  end subroutine printsummary
END MODULE PRINTING
  
MODULE ARGUMENTS
  use param
CONTAINS
  function getargval(arg) result(val)
    character(*), intent(in) :: arg
    integer :: val, n, p
    n = len_trim(arg)
    val = -1
    if (n == 0) return
    p = index('KMG', arg(n:n))
    if (p > 0) n = n-1
    read(arg(1:n),*) val
    if (2**(10*p) > huge(0)/val) stop 'Value must be < 2G'
    val = val*2**(10*p)
  end function getargval
  
  subroutine getarglist(arg, list)
    character(*), intent(in) :: arg
    integer, allocatable, intent(inout) :: list(:)
    integer :: j=1, n
    if (len_trim(arg) == 0) return
    n = min(19, count([(arg(j:j),j=1,len_trim(arg))] == ',')) + 1
    allocate(list(n))
    read(arg,*) list
  end subroutine getarglist
  
  logical function getargs(s, memsiz, n, nthr)
    ! Return true if there are any arguments
    real(dble), intent(out)              :: s      ! seconds
    integer,    intent(out)              :: memsiz ! memory size
    integer,    intent(out), allocatable :: n(:)   ! matrix sizes
    integer,    intent(out), allocatable :: nthr(:)! number of threads
    !
    integer :: nargs, i
    character :: arg*40
    !
    nargs = command_argument_count()
    getargs = (nargs > 0)
    if (.not.getargs) return
    !
    s = -1
    memsiz = -1
    do i = 1,nargs
      call get_command_argument(i, arg)
      select case(arg(1:min(2,len(arg))))
      case('-n')
        call getarglist(arg(3:), n)
      case('-t')
        call getarglist(arg(3:), nthr)
      case('-s')
        s = getargval(arg(3:))
      case('-m')
        memsiz = getargval(arg(3:))
      case('-')
      case default
        print *, 'unkown option:', arg
        stop
      end select
    end do
    if (.not. allocated(n)) then
      n = [100, 1000]
    end if
    if (.not. allocated(nthr)) then
      nthr = [1]
    endif
    if (s < 0) s = 1
  end function getargs
END MODULE ARGUMENTS

PROGRAM DRIVER
  use timing
  use printing
  use arguments
  !
  implicit none
  !
  integer, parameter :: n_routines_to_test = 4
  !
  type(timing_info) :: time
  !
  character(5), parameter   :: blasrouts(6) = &
    & ['dgemv', 'dgemm', 'dsyrk', 'dtrsm', 'dtrmv', 'dgbmv']
  !                                           dgemv dgemm dsyrk dtrsm dtrmv dgbmv
  integer,      parameter   :: xdim(6)      = [  1,    2,    2,    2,    1,    1 ]     
  integer,      parameter   :: ydim(6)      = [  1,    2,    0,    2,    0,    1 ]
  integer,      parameter   :: blaslevel(6) = [  2,    3,    3,    3,    2,    2 ]
  real(dble),   parameter   :: fctr(6)      = [1d0,  1d0,  2d0,  2d0,  1d0, 1.8d0]
  integer,      allocatable :: nv(:), nthr(:)
  logical                   :: dotime
  integer                   :: memsiz, maxwskip, nw
  real(dble), allocatable   :: w(:)
  integer, allocatable      :: skips(:,:,:), wskip(:,:)
  real(dble)                :: seconds, allottedtime, factor, totaltime
  integer                   :: n, in, nn, nrout, irt, irout, nspread, j
  !
  external dgemv, dgemm, dsyrk, dtrsm, dtrmv, dgbmv
  !
  print *
  dotime = getargs(seconds, memsiz, nv, nthr)
  if (.not. dotime) then
    call printhelp()
  else
    totaltime = 0
    nn = ubound(nv,1)
    if (maxval(nv)*3 >= (2**30/maxval(nv))) then
      print '("Max matrix size is ", I0)', int(sqrt(2d0**31/6))
      stop
    end if
    nrout = ubound(blasrouts,1)
    allottedtime = seconds/12
    allocate(skips(nn, 3, nrout))
    allocate(wskip(nn, nrout))
    do irt = 1,nrout
      skips(:,:,irt) = reshape([nv**2, nv**xdim(irt), nv**ydim(irt)], [nn, 3])
    enddo
    wskip = 2*sum(skips,2)
    maxwskip = maxval(wskip)
    nw = max(memsiz, maxwskip)
    allocate (w(nw))
    call initw(w,nw)
    call print('header')
    ! Time once to minimize spinup effect
    call set_num_threads(nthr(1))
    n = nv(1)
    time = timeblas(blasrouts(1), w, nw, n, skips(1,:,1), wskip(1,1), allottedtime)
    do in = 1,nn ! e.g. loop n=100, 1000, 10000
      n = nv(in)
      do irout = 1,n_routines_to_test
        nspread = memsiz/wskip(in,irout)
        factor = 1d12/(real(n,dble)**blaslevel(irout)/fctr(irout))
        do j = 1,ubound(nthr,1)
          call set_num_threads(nthr(j))
          time = timeblas(blasrouts(irout), w, nw, n, skips(in,:,irout),&
            & wskip(in,irout), allottedtime)
          totaltime = timesum(time) + totaltime
          call print(blasrouts(irout), n, time, factor, nthr(j), nspread)
        enddo
      enddo
    enddo
    print '(A)', 'The table shows run time in picoseconds/BLAS-flop (flop is one * and one +)'
    call printsummary(memsiz, totaltime)
  endif
END PROGRAM DRIVER
