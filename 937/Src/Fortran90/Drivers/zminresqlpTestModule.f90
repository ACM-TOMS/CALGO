!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! File zminresqlpTestModule.f90
!
! This file illustrates how MINRESQLP can call Aprod, Aprodmtx, or Msolve
! with a short fixed parameter list, even if it needs arbitrary other data.
!
!
! Authors:
!     Sou-Cheng Choi <sctchoi@uchicago.edu>
!     Computation Institute (CI)
!     University of Chicago
!     Chicago, IL 60637, USA
!
!     Michael Saunders <saunders@stanford.edu>
!     Systems Optimization Laboratory (SOL)
!     Stanford University
!     Stanford, CA 94305-4026, USA
!
! History:
! 12 Jul 2011: Created complex version zminresqlpTestModule.f90
!              from real version minresqlpTestModule.f90.
! 16 Sep 2012: zAprodmtx written for Hermitian A in CCH format from Matrix Market.
! 28 Oct 2012: Debugged mm_ioModule.f90.
!              mm_file_read must read complex values as two doubles.
!              We are now able to test "CCH" matrices from Matrix Market
!              (representation type = coordinate, complex, Hermitian).
! 30 Oct 2012: shift declared real(dp), not complex(dp).
!              zxcheck in zminresqlpCheckModule now used to check x from zMINRESQLP.
! 05 Jan 2013: Use minresqlpReadMtxModule for both real and complex MM mtx files.
!              Matrix Market arrays are now allocatable.
!              AprodMtxCCH replaces zAprodmtx (just the name).
! 09 Sep 2013: Removed unused references to dnrm2, present, precon, relTol
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module zminresqlpTestModule

  use  zminresqlpDataModule,    only : zzero, zone, dp, sp, ip, prcsn, eps, zero, debug
  use  zminresqlpModule,        only : ZMINRESQLP, ZSYMORTHO
  use  zminresqlpBlasModule,    only : znrm2
  use   minresqlpReadMtxModule, only : ReadMtxSize, ReadMtx, nnzmax

  implicit none

  private                                       ! sets default for module
  public   :: zminresqlptest, zminresqlptest2, zminresqlptestMtxCCH, zsymorthotest
  private  :: Aprod, Msolve, zAprod, zMsolve, AprodMtxCCH, zxcheck

  ! DYNAMIC WORKSPACE DEFINED HERE.
  ! It is allocated in zminresqlptest and used by zAprod or zMsolve.

  real(dp),    allocatable :: d(:)     ! Defines diagonal matrix D for Aprod.
  complex(dp), allocatable :: d1(:)    ! Defines vector zAprod.
  real(dp)                 :: Ashift   ! Shift diagonal elements of D in  Msolve.
  real(dp)                 :: Mpert    ! Perturbation to D in Msolve
                                                ! to avoid having an exact preconditioner.

  integer(ip)              :: nnz      ! These ones are used by the
  integer(ip), allocatable :: indx(:)  ! Matrix Market Aprod routine
  integer(ip), allocatable :: jndx(:)  ! AprodMtxCCH
  real(dp)   , allocatable :: dval(:)
  integer(ip), allocatable :: ival(:)
  real(sp)   , allocatable :: rval(:)
  complex(dp), allocatable :: cval(:)

contains

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine Aprod (n,x,y)

    integer(ip), intent(in)    :: n
    complex(dp), intent(in)    :: x(n)
    complex(dp), intent(out)   :: y(n)

    !-------------------------------------------------------------------
    ! Aprod  computes y = A*x for some matrix A.
    ! This is a simple example for testing zMINRESQLP.
    ! A = diag(d), where d is a real vector.
    !-------------------------------------------------------------------

    integer(ip) :: i

    do i = 1, n
       y(i) = d(i)*x(i)
    end do

  end subroutine Aprod

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine zAprod (n,x,y)

    integer(ip), intent(in)    :: n
    complex(dp), intent(in)    :: x(n)
    complex(dp), intent(out)   :: y(n)

    !-------------------------------------------------------------------
    ! zAprod  computes y = A*x for some matrix A.
    ! This is a simple example for testing MINRESQLP.
    ! A has complex vectors d1 and conjg(d1) on its first super- and sub-diagonal.
    !-------------------------------------------------------------------

    intrinsic    :: conjg, size
    integer(ip)  :: i

    if (size(x) > 1) then
       y(1) = d1(1)*x(2)
    end if

    do i = 2, n-1
       y(i) = conjg(d1(i-1))*x(i-1) + d1(i)*x(i+1)
    end do

    if (size(x) > 1) then
       y(n) = conjg(d1(n-1))*x(n-1)
    end if

  end subroutine zAprod

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine Msolve(n,x,y)

    integer(ip), intent(in)    :: n
    complex(dp), intent(in)    :: x(n)
    complex(dp), intent(out)   :: y(n)

    !-------------------------------------------------------------------
    ! Msolve solves M*y = x for some Hermitian positive-definite matrix M.
    ! This is a simple example for testing MINRESQLP.
    ! Ashift will be the same as shift in MINRESQLP.
    !
    ! If Mpert = 0, the preconditioner will be exact, so
    ! zMINRESQLP should require either one or two iterations,
    ! depending on whether (A - shift*I) is positive definite or not.
    !
    ! If Mpert is nonzero, somewhat more iterations will be required.
    !-------------------------------------------------------------------

    intrinsic   :: abs, mod
    integer(ip) :: i
    real(dp)    :: di

    do i = 1, n
       di   = d(i) - Ashift
       if (debug) write(*,*) "i = ", i, ", di = ", di, ", di = ", d(i), ", Ashift = ", Ashift
       if (mod(i,10) == 0) di = di + Mpert
       if (abs(di) > eps) then
          y(i) = x(i) / di
       else
          y(i) = x(i)
       end if
    end do

  end subroutine Msolve

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine zMsolve(n,x,y)

    integer(ip), intent(in)    :: n
    complex(dp), intent(in)    :: x(n)
    complex(dp), intent(out)   :: y(n)

    !-------------------------------------------------------------------
    ! zMsolve solves M*y = x for some Hermitian positive-definite matrix M.
    ! This is a simple example for testing MINRESQLP.
    ! Ashift will be the same as shift in MINRESQLP.
    !-------------------------------------------------------------------

    y = x   ! M is identity for now

  end subroutine zMsolve

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine AprodMtxCCH (n,x,y)

    integer(ip), intent(in)  :: n
    complex(dp), intent(in)  :: x(n)
    complex(dp), intent(out) :: y(n)

    !-------------------------------------------------------------------
    ! AprodMtxCCH computes y = A*x for some Hermitian matrix A.
    ! stored in Matrix Market CCH format (coordinate complex hermitian).
    ! Only subdiagonal and diagonal elements are in (indx, jndx, cval).
    !
    ! 16 Sep 2012: Diagonals treated as REAL here to guard against
    !              CCS format (coordinate complex symmetric), which
    !              might have complex diagonals.
    !-------------------------------------------------------------------

    intrinsic    :: conjg
    integer(ip)  :: i, j, k
    complex(dp)  :: d

    y(1:n) = zzero

    do k = 1, nnz
       i = indx(k)
       j = jndx(k)
       d = cval(k)

       if (i > j) then          ! d = subdiagonal
          y(i) = y(i) + d*x(j)
          y(j) = y(j) + conjg(d)*x(i)
       else                     ! i = j, d = diagonal
          y(i) = y(i) + real(d)*x(i)
       end if

       ! if (k <= 10) write(*,*) '  ', i, ' ', j, ' ', d,  ' ', y(i)
    end do

  end subroutine AprodMtxCCH

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine zminresqlptest( n, precon, shift, pertM, sing, consis, default, nout )

    integer(ip), intent(in)    :: n, nout
    logical,     intent(in)    :: precon, sing, consis
    real(dp),    intent(in)    :: shift, pertM
    logical,     intent(in), optional :: default

    !-------------------------------------------------------------------
    ! minresqlptest solves sets up and solves a system (A - shift*I)x = b,
    ! using Aprod to define A and Msolve to define a preconditioner.
    !-------------------------------------------------------------------

    intrinsic :: real

    complex(dp) :: b(n), r1(n), w(n), x(n), xtrue(n), y(n)
    logical     :: disable, use_default
    integer(ip) :: j, itnlim, istop, itn, nprint
    real(dp)    :: Anorm, Acond, Arnorm, rnorm, rtol, r1norm, xnorm
    real(dp)    :: enorm, etol, wnorm, xnormtrue
    real(dp)    :: maxxnorm, TranCond, Acondlim

    character(len=*), parameter ::  headerStr =                                   &
       "(// '-----------------------------------------------------'"    //        &
       "  / 'Test of zMINRESQLP.'"                                      //        &
       "  / '-----------------------------------------------------'"    //        &
       "  / 'shift  =', f12.4, 6x, 'pertM  =', f12.4)"
    character(len=*), parameter ::  footerStr1 =                                  &
       "(// ' zminresqlp appears to be successful.  n =', i7, '  Itns =', i7," // &
       "'  Relative error in x =', 1p, e8.1)"
    character(len=*), parameter ::  footerStr2 =                                  &
       "(// ' zminresqlp appears to have failed.    n =', i7, '  Itns =', i7," // &
       "'  Relative error in x =', 1p, e8.1)"
    character(len=*), parameter ::  debugStr1 =                                   &
       "(// 'Final residual =', 1p, e8.1)"
    character(len=*), parameter ::  debugStr2 =                                   &
       "(// 'Solution  x', 1p, 4(/e14.6, ' + I * ', e14.6))"

    write(nout, headerStr) shift, pertM

    allocate( d(n)  )          ! real array used in Aprod and Msolve to define A and M.
    allocate( d1(n) )          ! complex array used in zAprod to define A.
    Ashift = shift
    Mpert  = pertM

    if (.not. sing) then
      do j = 1, n                     ! Set d(*) to define A.
         d(j) = real(j,dp)/real(n,dp) ! We don't want exact integers.
         xtrue(j) = real(n+1-j,dp)    ! Set the true solution and the rhs
      end do                          ! so that (A - shift*I)*xtrue = b.
    else
      do j = 1, n-2                   ! Set d(*) to define A.
         d(j) = real(j,dp)/real(n,dp) ! We don't want exact integers.
         xtrue(j) = real(n+1-j,dp)    ! Set the true solution and the rhs
      end do                          ! so that (A - shift*I)*xtrue = b.
      d(n-1:n) = zero
      xtrue(n-1:n) = zzero
    end if

    call Aprod (n,xtrue,b)   ! b = A*xtrue
    b      = b - shift*xtrue ! Now b = (A - shift*I)*xtrue
    if (.not. consis) then
       b(n-1:n) = zone       ! b not in the range of A
    end if

    !debug = .false.
    if (debug) then
       write(*,*)
       write(*,*) 'A = diag(d), d = '
       do j = 1, n
          write(*,*) d(j)
       end do
       write(*,*) 'b = '
       do j = 1, n
          write(*,*) b(j)
       end do
    end if

    ! Set other parameters and solve.
    disable  = .false.
    itnlim   = n*3
    rtol     = 1.0e-12_dp
    maxxnorm = 1.0e+6_dp
    TranCond = 1.0e+7_dp
    Acondlim = 1.0e+15_dp

    if (prcsn == 6) then                    ! single precision is used
        rtol     = 1e-6_dp                  ! use a bigger rtol
        maxxnorm = 1.0e+4_dp                ! use a smaller maxxnorm
        TranCond = 1.0e+3_dp
        Acondlim = 1.0e+7_dp
    end if


    if (debug) then
       write(*,*)
       write(*,*)  'n = ', n, ', precon = ', precon, ', shift = ', shift, ', pertM = ', pertM, &
                   ', sing = ', sing, ', consis = ', consis, ', nout = ', nout
       write(*,*)
       write(*,*)  'itnlim = ', itnlim, ', nout = ', nout,                                     &
                   ', maxxnorm = ', maxxnorm, ', TranCond = ', TranCond, ', Acondlim = ', Acondlim
    end if

    if (present(default)) then
        use_default = default
    else
        use_default = .false.
    end if

    if (use_default) then
       call ZMINRESQLP( n=n, Aprod=Aprod, b=b, shift=shift, nout=nout, x=x, itn=itn )
    else
      if (precon) then
         call ZMINRESQLP( n, Aprod, b, shift, Msolve,  disable,                 &
                          nout, itnlim, rtol, maxxnorm, trancond, Acondlim,     &
                          x, istop, itn, rnorm, Arnorm, xnorm, Anorm, Acond )
      else
         call ZMINRESQLP( n=n, Aprod=Aprod, b=b, shift=shift, disable=disable,    &
                          nout=nout, itnlim=itnlim, rtol=rtol, maxxnorm=maxxnorm, &
                          trancond=trancond, Acondlim=Acondlim, x=x, istop=istop, &
                          itn=itn, rnorm=rnorm, Arnorm=Arnorm, xnorm=xnorm,       &
                          Anorm=Anorm, Acond=Acond )
       end if
    end if

    if (debug) then
       call Aprod (n,x,y)       ! y = A*x
       r1     = b - y + shift*x ! Final residual r1 = b - (A - shift*I)*x.
       r1norm = znrm2(n,r1,1)
       write(nout,debugStr1) r1norm

       nprint = min(n,20)
       write(nout,*) ' '
       write(nout,*) 'Some of x'
       write(nout,debugStr2) (j, real(x(j),dp), aimag(x(j)), j=1,nprint)
    else
       nprint = min(n,5)
       if (debug) then
         write(nout,*) ' '
         write(nout,*) 'Some of b and x'
         do j=1,nprint
            write(nout,*) b(j), x(j)
         end do
       end if
    end if

    w         = x - xtrue         ! Print a clue about whether
    wnorm     = znrm2(n,w,1)      ! the solution looks OK.
    xnormtrue = znrm2(n,xtrue,1)
    enorm     = wnorm/xnormtrue
    etol      = 1.0e-5_dp

    if (prcsn == 6) then                    ! single precision is used
          etol = 1e-3_dp                    ! use a bigger etol
    end if
    
    if (enorm <= etol) then
       write(nout, footerStr1) n, itn, enorm
    else
       write(nout, footerStr2) n, itn, enorm
    end if

    deallocate(d,d1)              ! Free work arrays

  end subroutine zminresqlptest

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine zminresqlptest2( n, precon, shift, pertM, nout )

    integer(ip), intent(in) :: n, nout
    logical,     intent(in) :: precon
    real(dp),    intent(in) :: shift, pertM

    !-------------------------------------------------------------------
    ! zminresqlptest2 solves sets up and solves a system (A - shift*I)x = b,
    ! using zAprod to define A and zMsolve to define a preconditioner.
    !-------------------------------------------------------------------

    intrinsic :: real

    complex(dp)  :: b(n), r1(n), w(n), x(n), xtrue(n), y(n)
    logical      :: disable
    integer(ip)  :: j, itnlim, istop, itn, nprint
    real(dp)     :: Anorm, Acond, Arnorm, rnorm, rtol, r1norm, xnorm
    real(dp)     :: enorm, etol, wnorm, xnormtrue
    real(dp)     :: maxxnorm, TranCond, Acondlim

    character(len=*), parameter ::  headerStr =                                   &
       "(// '-----------------------------------------------------'"    //        &
       "  / 'Test of zMINRESQLP.'"                                      //        &
       "  / '-----------------------------------------------------'"    //        &
       "  / 'shift  =', f12.4, 6x, 'pertM  =', f12.4)"
    character(len=*), parameter ::  footerStr1 =                                  &
       "(// ' zminresqlp appears to be successful.  n =', i7, '  Itns =', i7," // &
       "'  Relative error in x =', 1p, e8.1)"
    character(len=*), parameter ::  footerStr2 =                                  &
       "(// ' zminresqlp appears to have failed.    n =', i7, '  Itns =', i7," // &
       "'  Relative error in x =', 1p, e8.1)"
    character(len=*), parameter ::  debugStr1 =                                   &
       "(// 'Final residual =', 1p, e8.1)"
    character(len=*), parameter ::  debugStr2 =                                   &
       "(// 'Solution  x', 1p, 4(i6, e14.6, ' + I * ', e14.6))"

    write(nout, headerStr) shift, pertM

    allocate( d1(n-1) )          ! Array used in zAprod and zMsolve to define A and M.
    Ashift = shift
    Mpert  = pertM

    do j = 1, n-1               ! Set d(*) to define A.
       d1(j) = cmplx(real(j,dp)/real(n,dp),real(j,dp)/real(n,dp),dp)
                                ! We don't want exact integers.
    end do
    xtrue(1) = cmplx(0.550561797752808_dp, 2.707865168539327_dp, dp)
    xtrue(2) = cmplx(1.999999999999999_dp, 2.000000000000004_dp, dp)
    xtrue(3) = cmplx(2.146067415730337_dp, 2.775280898876405_dp, dp)
    xtrue(4) = cmplx(4.000000000000004_dp, 4.000000000000002_dp, dp)
    xtrue(5) = cmplx(5.168539325842692_dp, 4.359550561797752_dp, dp)

    call zAprod (n,xtrue,b)     ! b = A*xtrue
    b      = b - shift*xtrue    ! Now b = (A - shift*I)*xtrue

    !debug = .true.
    if (debug) then
       write(*,*)
       write(*,*) 'A = diag(d1,1) + conjg(diag(d1,-1)), d1 = '
       do j = 1, n-1
          write(*,*) d1(j)
       end do
       write(*,*) 'b = '
       do j = 1, n
          write(*,*) b(j)
       end do
    end if

    ! Set other parameters and solve.
    disable  = .false.
    itnlim   = n*3
    rtol     = 1.0e-12_dp
    maxxnorm = 1.0e+7_dp
    TranCond = 1.0e+7_dp
    Acondlim = 1.0e+15_dp

    if (prcsn == 6) then                    ! single precision is used
        rtol     = 1e-4_dp                  ! use a bigger rtol
        maxxnorm = 1.0e+4_dp                ! use a smaller maxxnorm
        TranCond = 1.0e+3_dp
        Acondlim = 1.0e+7_dp
    end if

    if (debug) then
       write(*,*)
       write(*,*) 'n = ', n, ', precon = ', precon, ', shift = ', shift, ', pertM = ', pertM, &
                  ', nout = ', nout
       write(*,*)
       write(*,*) 'itnlim = ', itnlim, ', nout = ', nout,                                     &
                  ', maxxnorm = ', maxxnorm, ', TranCond = ', TranCond, ', Acondlim = ', Acondlim
    end if

    call ZMINRESQLP( n, zAprod, b, shift, zMsolve, disable,               &
                     nout, itnlim, rtol, maxxnorm, trancond, Acondlim,    &
                     x, istop, itn, rnorm, Arnorm, xnorm, Anorm, Acond )
    call zAprod (n,x,y)       ! y = A*x
    r1     = b - y + shift*x  ! Final residual r1 = b - (A - shift*I)*x.
    r1norm = znrm2(n,r1,1)
    if (debug) write(nout,debugStr1) r1norm

    nprint = min(n,20)
    if (debug) &
      write(nout,debugStr2) (j, real(x(j),dp), aimag(x(j)), j=1,nprint)  ! Print some of the solution

    w      = x - xtrue                      ! Print a clue about whether
    wnorm  = znrm2(n,w,1)                   ! the solution looks OK.
    xnormtrue = znrm2(n,xtrue,1)
    enorm  = wnorm/xnormtrue
    etol   = 1.0e-5_dp
    if (enorm <= etol) then
       write(nout, footerStr1) n, itn, enorm
    else
       write(nout, footerStr2) n, itn, enorm
    end if

    deallocate(d1)                          ! Free work array

  end subroutine zminresqlptest2

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine zminresqlptestMtxCCH(input_file, consis, nout, tol)

    character(80),  intent(in)           :: input_file
    logical,        intent(in)           :: consis
    integer(ip),    intent(in)           :: nout
    real(dp),       intent(in), optional :: tol

    ! 30 Oct 2012: Use zxcheck to check computed x from zminresqlp.
    ! 02 Jan 2013: Print n in the "zminresqlp  appears to be ..." message
    !              to help identify the problem.

    intrinsic      :: real

    integer(ip)    :: input_unit, nrow, ncol
    character(14)  :: id
    character(10)  :: rep
    character( 6)  :: type
    character(7)   :: field
    character(19)  :: symm

    complex(dp), allocatable  :: b(:), x(:), r1(:), w(:)

    logical     :: disable
    integer(ip) :: n, j, itn, itnlim, istop, inform
    real(dp)    :: shift, Anorm, Acond, Arnorm, rnorm, rtol, xnorm
    real(dp)    :: maxxnorm, TranCond, Acondlim
    real(dp)    :: realj, realn, test1, test2, tolcheck

    character(len=*), parameter :: headerStr =              &
       "(// '----------------------------------------'"  // &
       "  / ' Test of zMINRESQLP on an MM CCH matrix '"  // &
       "  / '----------------------------------------')"

    character(len=*), parameter ::  footerStr1 = &
       "(// ' zminresqlp appears to be successful.  n =', i7, '  Itns =', i7," // &
       "'  test(r) =', 1p, e9.2, '  test(Ar) =', 1p, e9.2)"
    character(len=*), parameter ::  footerStr2 = &
       "(// ' zminresqlp appears to have failed.    n =', i7, '  Itns =', i7," // &
       "'  test(r) =', 1p, e9.2, '  test(Ar) =', 1p, e9.2)"
!    character(len=*), parameter ::  debugStr1 = &
!       "(// 'Final residual =', 1p, e8.1)"
!    character(len=*), parameter ::  debugStr2 = &
!       "(// 'Solution  x', 4(i6, e14.6))"

    write(nout, headerStr)

    call ReadMtxSize( input_file, input_unit, &
                      id, type, rep, field, symm, nrow, ncol, nnz )
    rewind( input_unit )

    ! Now we know the size of the problem.
    ! We should allocate only the arrays that will be used by the MM routines.
    ! For simplicity we allocate them all.

    nnzmax = nnz
    allocate( indx(nnz), jndx(nnz), ival(nnz) )
    allocate( rval(nnz) )   ! CRS
    allocate( dval(nnz) )   ! CDS
    allocate( cval(nnz) )   ! CCH

    call ReadMtx( input_file, input_unit, &
                  id, rep, field, symm, nrow, ncol, nnz, &
                  indx, jndx, ival, rval, dval, cval )

    n = nrow

    allocate( b(n) )
    allocate( x(n) )
    allocate( r1(n) )
    allocate( w(n) )

    realn = real(n,dp)
    do j = 1, n
       realj = real(j,dp)
       x(j)  = cmplx(realj,realj,dp) / realn
       b(j)  = zone
    end do

    if (debug) then
      write(nout,*) 'consis     = ', consis
      write(nout,*) 'input_file = ', trim(input_file)
      write(nout,*) 'n = ',  n, '  nnz = ', nnz
    end if

    ! Set other parameters and solve.
    disable  = .false.
    itnlim   = n*50
    rtol     = 1.0e-7_dp
    maxxnorm = 1.0e+8_dp
    TranCond = 1.0e+8_dp
    Acondlim = 1.0e+15_dp
    shift    = zero

    if (prcsn == 6) then                    ! single precision is used
        rtol     = 1e-4_dp                  ! use a bigger rtol
        maxxnorm = 1.0e+4_dp                ! use a smaller maxxnorm
        TranCond = 1.0e+3_dp
        Acondlim = 1.0e+7_dp
    end if

    if (consis) then
      call AprodMtxCCH (n,x,b)   ! b = A*x
      if (debug) then
        write(nout,*) ' '
        write(nout,*) 'norm(b) =', znrm2(n,b,1)
        write(nout,*) 'Some of the x defining b'
        do j = 1, min(n,5)
           write(nout,*) j, x(j)
        end do
      end if
    end if

    if (debug) then
      write(nout,*) 'Some of b'
      do j = 1, min(n,5)
         write(nout,*) j, b(j)
      end do
    end if

    if (debug) then
       write(*,*)
       write(*,*)  'n = ', n, ', shift = ', shift,  ', consis = ', consis, ', nout = ', nout
       write(*,*)
       write(*,*)  'itnlim = ', itnlim, ', nout = ', nout,                                       &
                   ', maxxnorm = ', maxxnorm, ', TranCond = ', TranCond, ', Acondlim = ', Acondlim
    end if

    call ZMINRESQLP( n=n, Aprod=AprodMtxCCH, b=b, shift=shift, disable=disable,&
                     nout=nout, itnlim=itnlim, rtol=rtol, maxxnorm=maxxnorm,   &
                     trancond=trancond, Acondlim=Acondlim, x=x, istop=istop,   &
                     itn=itn, rnorm=rnorm, Arnorm=Arnorm, xnorm=xnorm,         &
                     Anorm=Anorm, Acond=Acond )

    tolcheck = 10_dp*tol
    call zxcheck( n, AprodMtxCCH, b, shift, x, Anorm, tolcheck, nout, &
                  test1, test2, inform )

    if (inform <= 2) then
       write(nout, footerStr1) n, itn, test1, test2
    else
       write(nout, footerStr2) n, itn, test1, test2
    end if

    if (debug) then
      write(nout,*) ' '
      write(nout,*) 'norm(x) =', znrm2(n,x,1)
      write(nout,*) 'Some of the computed x'
      do j = 1, min(n,5)
         write(nout,*) x(j)
      end do
    end if

    deallocate( indx, jndx, ival, rval, dval, cval )
    deallocate( b, x, r1, w )

  end subroutine zminresqlptestMtxCCH

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine zsymorthotest( a, b, c_true, s_true, r_true, nout, tol)
    complex(dp), intent(in)           :: a, b, s_true, r_true
    real(dp),    intent(in)           :: c_true
    real(dp),    intent(in), optional :: tol
    integer(ip), intent(in)           :: nout

    !-------------------------------------------------------------------
    ! 20 Aug 2012: First version of zsymorthotest.
    !-------------------------------------------------------------------

    real(dp)                     :: c, norm_diff, relTol
    complex(dp)                  :: s, r, out(3), expected_out(3)

    character(len=*), parameter  :: footerStr1 = &
       "(// ' zsymortho  appears to be successful.  Relative error in [c,s,r] =', 1p, e8.1)"
    character(len=*), parameter  :: footerStr2 = &
       "(// ' zsymortho  appears to have failed.    Relative error in [c,s,r] =', 1p, e8.1)"

    call zsymortho( a, b, c, s, r )

    out          = (/ cmplx(     c,zero,dp), s     , r      /)
    expected_out = (/ cmplx(c_true,zero,dp), s_true, r_true /)
    norm_diff  = znrm2( 3, out - expected_out, 1 ) / znrm2( 3, expected_out, 1 )

    if (.not. present(tol)) then
       relTol = eps
    else
       relTol = tol
    end if

    write(nout,*)  ' '
    write(nout,*)  '-----------------------------------------------------'
    write(nout,*)  'Test of  zSYMORTHO.'
    write(nout,*)  '-----------------------------------------------------'
    write(nout,*)  'a = ',  a, '  b = ', b,  '  relTol = ', relTol
    write(nout,*)  ' '
    write(nout,*)  '[c,s,r]      = ',  out
    write(nout,*)  'true [c,s,r] = ',  expected_out

    if (norm_diff < relTol) then
       write(nout,footerStr1) norm_diff
    else
       write(nout,footerStr2) norm_diff
    end if

  end subroutine zsymorthotest

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! subroutine zxcheck
!
! zxcheck tests if a given x seems to be a solution of symmetric system
! Ax = b or Ax ~= b.
!
! 29 Oct 2012: derived from lsmrCheckModule.f90.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine zxcheck( n, Aprod, b, shift, x, Anorm, tol, nout, &
                      test1, test2, inform )

    integer(ip), intent(in)    :: n        ! No. of rows and cols of A
    integer(ip), intent(in)    :: nout     ! Output file number
    integer(ip), intent(out)   :: inform   ! = 0 if b = 0 and x = 0.
                                           ! = 1 or 2 if x seems to
                                           !   solve systems 1 or 2 below;
                                           ! = 4 otherwise.
    real(dp),    intent(in)    :: Anorm    ! An estimate of norm(A - shift*I)
                                           ! provided by zMINRESQLP.
    real(dp),    intent(in)    :: tol      ! tolerance for judging residuals.
                                           ! Typically the tol used for computing x.
    real(dp),    intent(in)    :: shift    ! Often zero.
    complex(dp), intent(in)    :: b(n)     ! The right-hand side of (A - shift*I)x ~= b.
    complex(dp), intent(in)    :: x(n)     ! The given solution estimate.
    real(dp),    intent(out)   :: test1    ! Should be small for consistent systems.
    real(dp),    intent(out)   :: test2    ! Should be small for LS systems.

    interface
       subroutine Aprod(n,x,y)          ! y := A*x
         use zminresqlpDataModule, only : dp, ip
         integer(ip),  intent(in)    :: n
         complex(dp),  intent(in)    :: x(n)
         complex(dp),  intent(out)   :: y(n)
       end subroutine Aprod
    end interface

    !-------------------------------------------------------------------
    ! One-liner: zxcheck tests if x solves Hermitian (A-shift*I)x=b.
    !
    ! Purpose:   zxcheck computes residuals and norms associated with
    ! the vector x and the problem solved by zMINRESQLP.
    ! It determines whether x seems to be a solution to either of the
    ! systems:  1. (A - shift*I)x = b
    !           2. min ||(A - shift*I)x - b||
    !
    ! History:
    ! 30 Oct 2012: zxcheck derived from xcheck in minresqlpCheckModule.f90.
    ! 02 Jan 2013: Print Anorm, xnorm, bnorm.
    !-------------------------------------------------------------------

    intrinsic           :: epsilon, max

    ! Local variables and arrays
    complex(dp)         :: r(n), v(n)
    real(dp)            :: bnorm, eps, rnorm, sigma, tol2, xnorm

    eps    = epsilon(eps)
    tol2   = max( tol, eps )

    call Aprod(n,x,r)        ! r = Ax
    r      = b - r + shift*x ! r = b - (A - shift*I)x = b - Ax + shift*x
    call Aprod(n,r,v)        ! v = Ar
    v      = v - shift*r     ! v = (A - shift*I)r = Ar - shift*r

    bnorm  = znrm2 (n,b,1)   ! Compute the norms of b, x, r, v.
    xnorm  = znrm2 (n,x,1)
    rnorm  = znrm2 (n,r,1)
    sigma  = znrm2 (n,v,1)

    if (debug) then
      if (nout > 0) write(nout,2200) shift, Anorm, xnorm, bnorm, rnorm, sigma
    end if

    !-------------------------------------------------------------------
    ! See if x seems to solve (A - shift*I)x = b  or  min ||(A - shift*I)x - b||.
    !-------------------------------------------------------------------
    if (bnorm == zero  .and.  xnorm == zero) then
       inform = 0
       test1  = zero
       test2  = zero
    else
       inform = 4
       test1  = rnorm / (Anorm*xnorm + bnorm)
       test2  = zero
       if (rnorm >  zero) test2  = sigma / (Anorm*rnorm)

       if (test2 <= tol2) inform = 2
       if (test1 <= tol2) inform = 1
    end if

    if (debug) then
      if (nout > 0) write(nout,3000) inform, tol2, test1, test2
    end if
    return

 2200 format(1p                    &
      // ' Enter zxcheck. Does x solve Ax = b?  where A is really (A - shift*I)' &
      /  '    shift     =', e10.3  &
      /  '    norm(A)   =', e10.3  &
      /  '    norm(x)   =', e10.3  &
      /  '    norm(b)   =', e10.3  &
      /  '    norm(r)   =', e15.8  &
      /  '    norm(Ar)  =', e10.3)
 3000 format(1p                    &
      /  '    inform    =', i2     &
      /  '    tol       =', e10.3  &
      /  '    test1     =', e10.3, ' (Ax = b)            rnorm/(Anorm*xnorm + bnorm)' &
      /  '    test2     =', e10.3, ' (least-squares)    Arnorm/(Anorm*rnorm)')

  end subroutine zxcheck

end module zminresqlpTestModule
