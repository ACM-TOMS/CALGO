program Decay

  ! Calculate amounts of radioactive isotopes for specified times, from given
  ! initial amounts

  use Connected_m, only: Connected
  use Dose_Factors_m, only: Dose_Factors
  use ENDF_m, only: ENDF
  use ISO_Fortran_Env, only: Error_Unit, Input_Unit, Output_Unit
  use Triangular_ODE, only: Compute_Solution, Compute_Z_Matrix
  use Triangularize_m, only: Triangularize

  implicit NONE

  integer, parameter :: RK = kind(1.0d0) ! Kind for real variables

  type :: What_t               ! Input data
    character(2) :: Sym = ''   ! Atomic symbol
    integer :: A = 0           ! Atomic weight
    integer :: Z = 0           ! Atomic number
    character :: M = ''        ! "m" or blank
    integer :: Row = 0         ! In activity matrix
    real(rk) :: AMU            ! Grams/mole
    real(rk) :: T_Half         ! Half life of decay mode, seconds
    real(rk) :: N0_g           ! Initial amount, grams
    real(rk) :: N0             ! Initial amount, moles
    real(rk) :: GBq = 0        ! Initial activity, Becquerels, from input
    character(2) :: D_Sym = '' ! Atomic symbol of daughter
    integer :: D_A = 0         ! Atomic weight of daughter
    integer :: D_Z = 0         ! Atomic number of daughter
    character :: D_M = ''      ! "m" or blank
    integer :: D_Row = 0       ! In activity matrix, of daughter
    real(rk) :: Branch = 0     ! Branch factor of decay mode
    real(rk) :: Dose = 0       ! Dose factor
  end type What_t

  type :: Scale_t              ! Scale factors for time
    character :: Code
    real(rk) :: Factor
  end type Scale_t

  real(rk), allocatable :: A(:,:)  ! Activity matrix, ln(2)/T_half
  real(rk), allocatable :: Amount(:,:)
  logical :: AnyDiff = .false.
  real(rk) :: Average           ! For checking mass conservation
  real(rk), parameter :: Avogadro = 6.022140857e23_rk
  real(rk), allocatable :: B(:,:)  ! Work area to permute A the lazy way
  integer, allocatable :: Blocks(:) ! Last row in each block of activity matrix
  integer :: C1, C2             ! Column numbers, for printing
  character(127) :: ERMSG
  logical :: Error = .false.
  real(rk), allocatable :: GBq_gm(:) ! Giga-Becquerels per gram for each row
  integer :: I
  integer :: In = Input_Unit
  integer :: IOSTAT
  character(127) :: IOMsg
  integer :: IScale = 1         ! Index in Scales, default seconds
  integer :: J                  ! Subscript to search ENDF
  integer :: K
  real(rk), parameter :: Log2 = log(2.0_rk)
  character(255) :: Line        ! of input
  logical :: LogTime = .false.  ! Step and Tn are base-10 log of time
  integer :: M
  logical :: Moles = .false.    ! Output moles instead of grams
  real(rk), allocatable :: More(:,:) ! for printing GBq/gm, dose factor
  integer :: N = 0              ! Line number of input
  integer :: NB = 0             ! Number of blocks of activity matrix
  integer :: NRow = 0           ! Number of rows (and columns) for A
  integer :: NStep              ! ceiling(t-finish - t-start) / t-step
  integer :: NWhat = 0          ! How much of What is used
  integer, allocatable :: NZ_Col(:) ! Nonzero columns of radiotoxicity
  integer :: Out = Output_Unit
  integer, allocatable :: P(:)  ! Permutation of A to make it triangular
  real(rk) :: Scale = 1.0       ! Scale for s m h d y times
  integer :: Stat
  real(rk) :: StdDev            ! For checking mass conservation
  real(rk) :: Step = 0.0        ! Time step or log10(time step)
  character(7) :: Sym, D_Sym    ! Input symbol and daughter symbol
  character(7), allocatable :: Syms(:) ! All the symbols
  real(rk) :: T_Half            ! For testing consistency with what(.)%GBq
  real(rk) :: T = 0.0           ! Time or log10(time)
  real(rk) :: Time
  real(rk) :: Tn = 0.0          ! Ending time or log10(ending time)
  integer :: Verbose = 0
  type(what_t), allocatable :: What(:), What_tmp(:)
  type(what_t) :: WR
  integer, allocatable :: Where(:), Where_tmp(:)  ! Index from row number to What
  integer :: W = 6              ! Width of activity array dump
  logical :: Zero = .false.     ! Print columns of zero radiotoxicity

  character(*), parameter :: InHead1 = &
    '       Isotope                                                '
  character(*), parameter :: InHead2 = &
    '                     Daughter          Branch    Dose'
  character(*), parameter :: InHead3 = &
    '       A   Sym   Z  Row   AMU        T_Half     N0 gm      '
  character(*), parameter :: InHead4 = &
    'N0 moles   GBq          A   Sym   Z  Row  Factor    Factor'
  character(*), parameter :: InFMT = & ! for printing input
    '(i3,":",i4,a,2x,a2,2i5,f11.5,1p,e12.4,3e11.3,i5,a,2x,a2,2i5,' // &
    '0p,f10.4,2x,1p,g7.1)'

  type(scale_t), parameter :: Scales(5) = [ &
    scale_t('s', 1.0), &
    scale_t('m', 60.0), &
    scale_t('h', 60.0 * 60.0), &
    scale_t('d', 60.0 * 60.0 * 24.0), &
    scale_t('y', 60.0 * 60.0 * 24.0 * 365.242190_rk ) ]

  allocate ( what(0:100), where(100) )

  ! Collect command-line arguments
  i = 0
  do
    i = i + 1
    call get_command_argument ( i, line )
    if ( line == '' ) exit
    select case ( line(1:2) )
    case ( '-l', '-L' )
      logTime = .true.
    case ( '-m', '-M' ) ! Moles
      moles = .true.
    case ( '-s', '-S' ) ! Starting time
      if ( line(3:) == '' ) then
        i = i + 1
        call get_command_argument ( i, line(3:) )
      end if
      read ( line(3:), *, iostat=stat, iomsg=iomsg ) step
      if ( stat /= 0 ) then
        write ( error_unit, '(a,i0,a/a/a)' ) 'Error ', stat, &
          & ' while trying to read step from', trim(line(3:)), trim(iomsg)
        call usage
      end if
    case ( '-t', '-T' ) ! Ending time
      if ( line(3:) == '' ) then
        i = i + 1
        call get_command_argument ( i, line(3:) )
      end if
      read ( line(3:), *, iostat=stat, iomsg=iomsg ) tn
      if ( stat /= 0 ) then
        write ( error_unit, '(a,i0,a/a/a)' ) 'Error ', stat, &
          & ' while trying to read ending time from', trim(line(3:)), trim(iomsg)
        call usage
      end if
    case ( '-u', '-U' ) ! Time units s m h d y
      if ( line(3:) == '' ) then
        i = i + 1
        call get_command_argument ( i, line(3:) )
      end if
      do iScale = 1, size(scales,1)
        if ( line(3:3) == scales(iScale)%code ) then
          scale = scales(iScale)%factor
          exit
        end if
      end do
      if ( iScale > size(scales,1) ) then
        write ( error_unit, '(a,99(1x,a))' ) 'Time units are not one of', scales%code
        call usage
      end if
    case ( '-v', '-V' )
      verbose = 1
      if ( line(3:) == '' ) then
        i = i + 1
        call get_command_argument ( i, line(3:) )
      end if
      if ( line(3:3) == '-' ) then
        i = i - 1
      else
        read ( line(3:), *, iostat=stat, iomsg=iomsg ) verbose
        if ( stat /= 0 ) then
          write ( error_unit, '(a,i0,a/a/a)' ) 'Error ', stat, &
            & ' while trying to read verbosity from', trim(line(3:)), trim(iomsg)
          call usage
        end if
      end if
    case ( '-w', '-W' )
      if ( line(3:) == '' ) then
        i = i + 1
        call get_command_argument ( i, line(3:) )
      end if
      read ( line(3:), *, iostat=stat, iomsg=iomsg ) w
      if ( stat /= 0 ) then
        write ( error_unit, '(a,i0,a/a/a)' ) 'Error ', stat, &
          & ' while trying to read width from', trim(line(3:)), trim(iomsg)
        call usage
      end if
    case ( '-z', '-Z' )
      zero = .true.
    case default
      if ( line(1:1) /= '-' ) exit
      call usage
    end select
  end do

  if ( tn == 0 ) then
    write ( error_unit, '(a)' ) 'No ending time specified'
    call usage
    stop
  end if

  if ( line /= '' ) then
    open ( 10, file=trim(line), status='old', form='formatted', &
         & access='sequential', iostat=iostat, iomsg=iomsg )
    if ( iostat /= 0 ) then
      write ( error_unit, '(a,i0,2a/a)' ) 'Error ', iostat, ' while opening ', &
        & trim(line), trim(iomsg)
      stop
    end if
    in = 10
  end if

  i = i + 1
  call get_command_argument ( i, line )
  if ( line /= '' ) then
    open ( 11, file=trim(line), form='formatted', &
         & access='sequential', iostat=iostat, iomsg=iomsg )
    if ( iostat /= 0 ) then
      write ( error_unit, '(a,i0,2a/a)' ) 'Error ', iostat, ' while opening ', &
        & trim(line), trim(iomsg)
      stop
    end if
    out = 11
  end if

  call get_command ( line )
  write ( out, '(2a)' ) 'Invoked using ',trim(line)

  if ( verbose > 1 ) then
    write ( out, '(/a)' ) 'Original input'
    write ( out, '(2a)' ) inhead1, inhead2
    write ( out, '(2a)' ) inhead3, inhead4
  end if

  do
    read ( in, '(a)', end=9 ) line
    n = n + 1
    if ( line == '' ) cycle
    j = scan(line,'!#%')
    if ( j /= 0 ) then
      if ( j == 1 .or. line(1:j-1) == '' ) cycle
    end if
    nWhat = nWhat + 1
    if ( nWhat > size(what) ) then ! Double size of What and Where
      allocate ( what_tmp(2*size(what)) )
      what_tmp(1:size(what)) = what
      call move_alloc ( what_tmp, what )
      allocate ( where_tmp(2*size(where)) )
      where_tmp(1:nRow) = where(1:nRow)
      call move_alloc ( where_tmp, where )
    end if
    ! Read a non-comment non-blank line
    read ( line, 1, iostat=iostat, iomsg=ermsg ) sym, d_sym, what(nWhat)%t_half, &
      & what(nWhat)%branch, what(nWhat)%N0_g, what(nWhat)%GBq
  1 format ( a7, t12, a7, t21, 4g10.0 )
    if ( iostat /= 0 ) then
      write ( error_unit, '(a,i0,a)' ) 'Error ', iostat, ' while reading', &
        & trim(line)
      nWhat = nWhat - 1
      error = .true.
      cycle
    end if
    if ( lookup ( sym, .false. ) ) then
      nWhat = nWhat - 1
      write ( error_unit, '(3a,i0,a)' ) 'Symbol ', trim(sym), ' at line ', n, &
        & ' not found'
      cycle
    end if
    if ( what(nWhat)%t_half /= 0 .and. what(nWhat)%d_sym /= 'stable' ) then
      if (  lookup ( d_sym, .true. ) ) then
        nWhat = nWhat - 1
        write ( error_unit, '(3a,i0,a)' ) 'Symbol ', trim(sym), ' at line ', n, &
          & ' not found'
        cycle
      end if
    else
      what(nWhat)%d_sym = what(nWhat)%sym
      what(nWhat)%d_a =  what(nWhat)%a
      what(nWhat)%d_z =  what(nWhat)%z
    end if
    what(nWhat)%d_a = abs(what(nWhat)%d_a)
    if ( verbose > 1 ) then
      write ( out, inFMT ) &
        & n, what(nWhat)%a, what(nWhat)%m, what(nWhat)%sym, what(nWhat)%z, nWhat, &
        & what(nWhat)%AMU, what(nWhat)%t_half, what(nWhat)%N0_g, what(nWhat)%N0, &
        & what(nWhat)%GBq, what(nWhat)%d_a, what(nWhat)%d_m, what(nWhat)%d_sym, what(nWhat)%d_z, &
        & what(nWhat)%d_row, what(nWhat)%branch, what(nWhat)%dose
    end if
  end do

9 continue

  ! Sort What increasing on A, then increasing on Z, then descending on M
  do i = 1, nWhat
    wr = what(i)
    j = i
    do while ( wr%a < what(j-1)%a .or. &
             & wr%a == what(j-1)%a .and. wr%z < what(j-1)%z .or. &
             & wr%a == what(j-1)%a .and. wr%z == what(j-1)%z &
               & .and. wr%m > what(j-1)%m )
      what(j) = what(j-1)
      j = j - 1
    end do
    what(j) = wr
  end do

  ! Calculate number of rows of activity matrix
  do i = 1, nWhat
    if ( what(i)%a /= what(i-1)%a .or. &
       & what(i)%z /= what(i-1)%z .or. &
       & what(i)%m /= what(i-1)%m ) then
      nRow = nRow + 1
      what(i)%row = nRow
      where(nRow) = i
    else
      what(i)%row = what(i-1)%row
    end if
  end do

  ! Calculate daughter row numbers
  do i = 1, nWhat
    do j = 1, nWhat
      if ( what(i)%d_sym == what(j)%sym .and. &
         & what(i)%d_a == what(j)%a .and. &
         & what(i)%d_m == what(j)%m ) exit
    end do
    if ( j <= nWhat ) then
      what(i)%d_row = what(j)%row
    else
      write ( error_unit, '(3a,i0,3a,i0,a,i0)' ) 'No row for daughter isotope ', &
        & trim(what(i)%d_sym), '-', what(i)%d_a, ' of ', &
        & trim(what(i)%sym), '-', what(i)%a, ' at row ', i
      error = .true.
    end if
  end do

  if ( error ) then
    write ( error_unit, '(a)' ) 'Errors in input'
    stop
  end if

  if ( verbose > 1 ) then
    write ( out, '(/a)' ) 'Sorted input'
    write ( out, '(2a)' ) inhead1, inhead2
    write ( out, '(2a)' ) inhead3, inhead4
    do i = 1, nWhat
      write ( out, inFMT )&
        & i, what(i)%a, what(i)%m, what(i)%sym, what(i)%z, what(i)%row, &
        & what(i)%AMU, what(i)%t_half, what(i)%N0_g, what(i)%N0, &
        & what(i)%GBq, what(i)%d_a, what(i)%d_m, what(i)%d_sym, what(i)%d_z, &
        & what(i)%d_row, what(i)%branch, what(i)%dose
    end do
  end if

  if ( verbose > 2 ) then
    ! Report differences between input and calculated half life
    do i = 1, nWhat
      if ( what(i)%GBq /= 0 ) then
        t_half = 1.0e-9 * log(2.0) * avogadro * what(i)%N0 / what(i)%GBq
        if ( abs(t_half - what(i)%t_half) > &
           & 0.5e-3 * ( t_half + what(i)%t_half ) ) then
          if ( .not. anyDiff ) write ( out, 6 )
          6 format (/'Half lives calculated from radioactivity are different from input' / &
                  & '               Calculated  Input       Relative' / &
                  & '     Isotope   T_half      T_Half      Difference')
          write ( out, '(i3,": ",a,"-",i0,a,t14,1p,3e12.4)' ) i, &
            & trim(what(i)%sym), what(i)%a, what(i)%m, t_half, &
            & what(i)%t_half, &
            & abs(t_half - what(i)%t_half) / ( 0.5 * ( t_half + what(i)%t_half ) )
          anyDiff = .true.
        end if
      end if
    end do
  end if

  ! Create the activity matrix
  allocate ( A(nRow,nRow), syms(nRow+1) )  
  a = 0.0
  do i = 1, nWhat
    if ( what(i)%t_half /= 0 ) then
      m = what(i)%row
      a(m,m) = a(m,m) - log2 / what(i)%t_half * what(i)%branch / 100
      j = what(i)%d_row
      a(j,m) = a(j,m) + log2 / what(i)%t_half * what(i)%branch / 100
    end if
  end do

  do j = 1, nRow
    write ( syms(j), '(a,"-",i0,a)' ) trim(what(where(j))%sym), &
          & what(where(j))%a, what(where(j))%m
  end do
  syms(nRow+1) = 'Total'

  if ( verbose > 2 ) then
    write ( out, '(/a)' ) 'Before permutation'
    call dumpMatrix ( 'Activity Matrix', a, [0,nRow] )
  end if

  ! Permute A to triangular form
  allocate ( p(nRow), blocks(0:nRow) )
  call connected ( a, p, blocks, nb )
  call triangularize ( a, p, blocks(0:nb) )
  allocate ( b(nRow,nRow) )
  do i = 1, nRow
    b(i:,i) = a(p(i:),p(i))
    b(i,i+1:) = a(p(i),p(i+1:))
  end do
  syms(1:nRow) = syms(p)

  ! Check that the permutation worked
  do i = 1, nRow
    if ( any(b(i,i+1:) /= 0) ) then
      write ( out, '(/ a,i0,3a,i0,2a)' ) 'Row ', i, ', ', trim(what(where(i))%sym), &
        & '-', what(where(i))%a, trim(what(where(i))%m), &
        & ', of activity matrix is not triangular'
      call dumpMatrix ( 'Activity Matrix', a, [0,nRow] )
      write ( out, '(a)' ) 'Block #   Permutation:'
      do k = 1, nb
        write ( out, '(i5,":",(t10,10i5))' ) k, p(blocks(k-1)+1:blocks(k))
      end do
      stop
    end if
  end do
  call move_alloc ( b, a )
  allocate ( b(nRow,nRow) )

  where(1:nRow) = where(p(1:nRow)) ! Permute the indices back to What too,
                                   ! so we can find stuff using row numbers

  ! Calculate radioactivity
  allocate ( GBq_gm(nRow) )
  do i = 1, nRow
    GBq_gm(i) = -1.0e-9_rk * avogadro * a(i,i) / log2 / what(where(i))%amu  
  end do

  if ( verbose > 1) then
    allocate ( more(nRow,2) )
    more(:,1) = GBq_gm
    more(:,2) = what(where(1:nRow))%dose
    if ( verbose > 2 ) &
      & write ( out, '(/a)' ) 'After permutation and partitioning into independent blocks'
    call dumpMatrix ( 'Activity Matrix', a, blocks(0:nb), &
      & ['GBq/gm','Dose  '], more )
  end if

  ! Check that column sums are zero
  do i = 1, nRow
    if ( sum(a(:,i)) > 2.0 * maxval(abs(a(:,i))) * epsilon(a(1,1)) ) then
      write ( error_unit, '(a,i0,a)' ) 'Sum of A(:,', i, &
        & ') is not zero.  Check branch factors.'
      stop
    end if
  end do

  ! Compute solution coefficients in B
  b = 0
  do i = 1, nb
    j = blocks(i-1) + 1
    k = blocks(i)
    call compute_Z_matrix ( a(j:k,j:k), b(j:k,j:k) )
  end do
  if ( verbose > 1 ) call dumpMatrix ( 'Solution Matrix', b(:,1:), blocks(0:nb) )

  ! Compute amounts
  nStep = ceiling( ( tn - t ) / step ) + 1
  allocate ( amount(nRow+1,nStep) )
  m = 0
  t = 0
  do while ( t <= tn )
    m = m + 1
    do i = 1, nb
      j = blocks(i-1) + 1
      k = blocks(i)
      time = merge(10**t,t,logTime) * scale
      if ( moles ) then
        call compute_solution ( time, a(j:k,j:k), b(j:k,j:k), &
          & what(where(j:k))%n0, amount(j:k,m) )
      else ! Don't compute with moles and multiply amount by
           ! single-precision AMU
        call compute_solution ( time, a(j:k,j:k), b(j:k,j:k), &
          & what(where(j:k))%n0_g, amount(j:k,m) )
      end if
    end do
    t = t + step
  end do

  ! Compute total amount
  do j = 1, m
    amount(nRow+1,j) = sum(amount(1:nRow,j))
  end do

  c2 = 0
  c1 = 1
  do while ( c1 <= nRow+1 )
    write ( out, '(/15x,a)' ) "Amounts, " // merge ( "moles","grams",moles)
    write ( out, '("Time ", a, 9x)', advance='no' ) scales(iScale)%code
    c2 = min(nRow+1, c1 + w - 1)
    call matrixHead ( 11, syms(c1:c2) )
    t = 0
    do j = 1, m
      time = merge(10**t,t,logTime) * scale
      write ( line, '(1pg13.5)' ) time / scales(iScale)%factor
      write ( out, '(a,t14,1p,99e11.3)' ) trim(adjustl(line)), amount(c1:c2,j)
      t = t + step
    end do
    c1 = c2 + 1
  end do

  average = sum(amount(nrow+1,1:m)) / m
  if ( verbose > 3 ) then
    ! Print total amount with 16 digit precision for checking mass conservation
    write ( out, '(/"Time ", a, t18, a, t40, a)' ) scales(iScale)%code, &
      & 'Total Amount', '( Amount - Average ) / Average'
    do j = 1, m
      time = merge(10**t,t,logTime) * scale
      write ( line, '(1pg13.5)' ) time / scales(iScale)%factor
      write ( out, '(a,t14,1p,2g25.18)' ) trim(adjustl(line)), amount(nRow+1,j), &
        & ( amount(nRow+1,j) - average ) / average
      t = t + step
    end do
  end if
  if ( m > 1 ) then
    stdDev = sqrt(sum((amount(nrow+1,1:m)-average)**2)/(m-1))
    write ( out, '(3(a,1x,1pe11.3))' ) 'Total Std Dev', stdDev, &
      & ' StdDev / Average ', stdDev / average
  end if

  ! Convert amounts to grams if printing moles
  if ( moles ) then
    do j = 1, m
      amount(1:nRow,j) = amount(1:nRow,j) * what(where(1:nRow))%amu
    end do
  end if

  ! Compute radiotoxicity in Sieverts
  do j = 1, m
    amount(1:nRow,j) = amount(1:nRow,j) * GBq_gm(1:nRow) * 1.0e9 * &
                     & what(where(1:nRow))%dose
    amount(nRow+1,j) = sum(amount(1:nRow,j))
  end do

  ! Determine which columns are nonzero
  allocate ( nz_col(nRow+1) )
  if ( zero ) then
    nz_col = [ ( k, k = 1, nRow+1 ) ]
    k = nRow + 1
  else
    k = 0
    do i = 1, nRow
      if ( any(amount(i,1:m) /= 0) ) then
        k = k + 1
        nz_col(k) = i
      end if
    end do
    k = k + 1
    nz_col(k) = nRow + 1
  end if

  ! Print radiotoxicity in Sieverts
  if ( k > 1 ) then
    c2 = 0
    c1 = 1
    do while ( c1 <= k )
      write ( out, '(/15x,"Radiotoxicity, Sieverts")' )
      write ( out, '("Time ", a, 9x)', advance='no' ) scales(iScale)%code
      c2 = min(k, c1 + w - 1)
      call matrixHead ( 11, syms(nz_col(c1:c2)) )
      t = 0
      do j = 1, m
        ! Print the time
        time = merge(10**t,t,logTime) * scale
        write ( line, '(1pg13.5)' ) time / scales(iScale)%factor
        write ( out, '(a,t14,1p,99e11.3)' ) trim(adjustl(line)), &
          & amount(nz_col(c1:c2),j)
        t = t + step
      end do
      c1 = c2 + 1
    end do
    if ( .not. zero ) &
      & write ( out, '(/a)' ) 'Columns with zero radiotoxicity not printed.'
  end if
  write ( out, '(/a)' ) 'Radiotoxicity is computed for adult ingestion using dose factors from ICRP'
  write ( out, '(a)' )  'publication 119, which does not provide dose factors for radioactive'
  write ( out, '(a)' )  'isotopes with half lives less than ten minutes or greater than 10^9 years.'
  write ( out, '(a)' )  'Dose factors for gases are given as Sv/day per Bq/m^3.  Radiotoxicity is'
  write ( out, '(a)' )  'not calculated for gases.'

contains

  subroutine DumpMatrix ( Title, A, Blocks, MoreLabel, More )
    character(len=*), intent(in) :: Title
    real(rk), intent(in) :: A(:,:)
    integer, intent(in) :: Blocks(0:)
    character(*), intent(in), optional :: MoreLabel(:)
    real(rk), intent(in), optional :: More(:,:) ! GBq/gm
    integer :: B       ! Index into Blocks
    integer :: C1, C2  ! Column/Row indices
    integer :: R, Top  ! Row subscript, top column subscript
    do b = 1, ubound(blocks,1)
      top = blocks(b)
      c2 = blocks(b-1)
      c1 = c2 + 1
      ! Print a block of A from row/column blocks(b-1)+1 to blocks(b)
      ! but with no more than W columns at a time.
      do while ( c1 <= top )
        c2 = min(top, c1 + w - 1)
        write ( out, '(/a)', advance='no' ) trim(title)
        write ( out, '(a)', advance='no' ) repeat(' ', 18-len_trim(title))
        call matrixHead ( 12, syms(c1:c2) )
        do r = blocks(b-1)+1, top
          write ( out, '(i3,"#",3x,a,t16,1p,99g12.4)' ) r, &
            & syms(r), a(r,c1:c2)
        end do
        if ( present(moreLabel) ) then
          do r = 1, size(moreLabel)
            write ( out, '(t8,a,t16,1p,99g12.4)' ) moreLabel(r), more(c1:c2,r)
          end do
        end if
        c1 = c2 + 1
      end do
    end do
  end subroutine DumpMatrix

  ! Find entries in dose factors, ENDF, and isotope tables
  logical function Lookup ( Sym, Daughter ) ! returns true if error
    character(7), intent(in) :: Sym
    logical, intent(in) :: Daughter
    integer :: A, I, J, K
    character :: M
    lookup = .false. ! Assume no error
    j = index(sym,'-')
    if ( j <= 1 ) then
      write ( error_unit, '(a,i0)' ) 'No hyphen in isotope identifier at line', &
        & n, trim(line)
      nWhat = nWhat - 1
      error = .true.
      lookup = .true.
      return
    end if
    k = scan(sym(j:),'mno')
    if ( k /= 0 ) then
      m = sym(j+k-1:j+k-1)
      if ( daughter ) then
        what(nWhat)%d_m = m
      else
        what(nWhat)%m = m
      end if
    else
      m = ''
      k = len_trim(sym)
    end if
    read ( sym(j+1:j+k-2), *, iostat=iostat, iomsg=ermsg ) a
    if ( iostat /= 0 ) then
      write ( out, '(a,i0,a)' ) 'Error ', iostat, &
        & ' while reading atomic weight from', trim(line)
      nWhat = nWhat - 1
      error = .true.
      lookup = .true.
      return
    end if
    if ( daughter ) then
      what(nWhat)%d_sym = sym(1:j-1)
      what(nWhat)%d_a = abs(a)
    else
      what(nWhat)%sym = sym(1:j-1)
      what(nWhat)%a = abs(a)
    end if
    do i = 1, size(ENDF)
      if ( ENDF(i)%sym == sym(1:j-1) .and. ENDF(i)%a == abs(a) .and. &
         & ENDF(i)%m == m ) exit
    end do
    if ( i > size(ENDF) ) then
      write ( out, '(2a)' ) 'No isotope information for ', trim(sym), trim(line)
      nWhat = nWhat - 1
      error = .true.
      lookup = .true.
      return
    end if
    if ( daughter ) then
      what(nWhat)%d_a = ENDF(i)%a
      what(nWhat)%d_z = ENDF(i)%z
    else
      what(nWhat)%a = ENDF(i)%a
      what(nWhat)%z = ENDF(i)%z
      ! Calculate number of moles
      what(nWhat)%amu = ENDF(i)%amu
      what(nWhat)%N0 = what(nWhat)%N0_g / what(nWhat)%amu
      ! Find dose factor
      do j = 1, size(dose_factors)
        if ( dose_factors(j)%sym == sym ) then
          what(nWhat)%dose = dose_factors(j)%dose_adult
          exit
        end if
      end do
    end if
  end function Lookup

  ! Column labels, same as the row labels
  subroutine MatrixHead ( Width, Syms )
    integer, intent(in) :: Width
    character(*), intent(in) :: Syms(:)
    integer :: I
    do i = 1, size(syms)-1
      line(1:width) = syms(i)
      write ( out, '(a)', advance='no' ) line(1:width)
    end do
    write ( out, '(a)' ) trim(syms(i))
  end subroutine MatrixHead

  subroutine Usage
    call get_command_argument ( 0, line )
    write ( out, '(3a)' )   'Usage: ', trim(line), ' [options] [input [output]]'
    write ( out, '(a)' )   ' Options:'
    write ( out, '(a)' )   '  -l     => Time step and Final time (-s and -t) are log10 of time.'
    write ( out, '(a)' )   '            If time step is not an integer, you might wish to'
    write ( out, '(a)' )   '            select a slightly larger final time, to compute the'
    write ( out, '(a)' )   '            desired number of steps.'
    write ( out, '(a)' )   '  -m     => Print moles instead of grams'
    write ( out, '(a)' )   '  -s[ ]# => Time step'
    write ( out, '(a)' )   '  -t[ ]# => Final time; always starts at zero'
    write ( out, '(a)' )   '  -u[ ]X => Time scale, X = s m h d y'
    write ( out, '(a)' )   '  -v[ ]# => Be verbose at level #'
    write ( out, '(a,i0)' ) '  -w[ ]# => Display # columns, default ', w
    write ( out, '(a)' )   '  -z     => Print columns having zero radiotoxicity'
    write ( out, '(a)' )   '  anything else => this output'
    write ( out, '(a)' )   'Option letters are not case sensitive.'
    write ( out, '(a)' )   'If input is not specified, input is from stdin'
    write ( out, '(a)' )   'If output is not specified, output is from stdout'
    stop
  end subroutine Usage

end program Decay
