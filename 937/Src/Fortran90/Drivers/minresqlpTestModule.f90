!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! File minresqlpTestModule.f90
!
! This file illustrates how MINRESQLP can call Aprod or Msolve
! with a short fixed parameter list, even if it needs arbitrary other data.
!
! Authors:
!     Sou-Cheng Choi <sctchoi@uchicago.edu>
!     Computation Institute
!     University of Chicago
!     Chicago, IL 60637, USA
!
!     Michael Saunders <saunders@stanford.edu>
!     Systems Optimization Laboratory (SOL)
!     Stanford University
!     Stanford, CA 94305-4026, USA
!
! History:
! 11 Oct 2007: First version of minresqlpTestModule.f90.
! 15 Oct 2007: Use real(8) everywhere.
! 16 Oct 2007: Use minresqlpDataModule to define dp = selected_real_kind(15).
! 12 Jul 2011: Created complex version zminresqlpTestModule.f90
!              from real version minresqlpTestModule.f90.
! 20 Aug 2012: Added Aprodmtx for multiplying a sparse matrix in
!              Matrix Market Format with a vector.
!              Also added two subroutines minresqlpmtxtest and symorthotest
!              for testing SYMORTHO and MINRESQLP on singular matrices
!              from the Matrix Market.
!              (Aprodmtx is now AprodMtxCPS, AprodMtxCRS, AprodMtxCDS.)
!              (minresqlpmtxtest is now
!               minresqlptestMtxCPS, minresqlptestMtxCRS, minresqlptestMtxCDS.)
! 08 Sep 2012: We are testing "CRS" matrices from Matrix Market
!              (representation type = coordinate, real, symmetric).
!              Found that the MM input routines for "real" type
!              store A in (indx, jndx, rval), not (indx, jndx, dval).
!              Simplified Aprodmtx now works.
! 29 Oct 2012: xcheck in minresqlpCheckModule now used to check x from MINRESQLP.
! 05 Jan 2013: Use minresqlpReadMtxModule for both real and complex MM mtx files.
!              Matrix Market arrays are now allocatable.
!              nnz is determined from the data file.
!              We don't have to guess nnzmax at compile time.
!              AprodMtxCRS replaces Aprodmtx (just the name).
!              AprodMtxCDS added to allow for MM matrices in CDS format (double).
! 21 Apr 2013: AprodMtxCPS added to allow for MM matrices in CPS format
!              (coordinate pattern symmetric, meaning nonzero Aij = 1.0).
!              minresqlptestMtxCRS replaces local subroutine minresqlpmtxtest.
!              minresqlptestMtxCDS, minresqlptestMtxCPS added.
! 09 Sep 2013: Removed unused references to realmin, dot_product, sqrt, 
!              present, precon, transpose, abs, sum.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module minresqlpTestModule

  use  minresqlpDataModule,    only : dp, sp, ip, prcsn, eps, zero, one, debug
  use  minresqlpModule,        only : MINRESQLP, SYMORTHO
  use  minresqlpBlasModule,    only : dnrm2
  use  minresqlpReadMtxModule, only : ReadMtxSize, ReadMtx, nnzmax

  implicit none

  private                                       ! sets default for module
  public   :: minresqlptest, &
              minresqlptestMtxCDS,  minresqlptestMtxCPS,  minresqlptestMtxCRS, &
              symorthotest
  private  :: Aprod, Msolve, &
              AprodMtxCDS, AprodMtxCPS, AprodMtxCRS, xcheck

  ! DYNAMIC WORKSPACE DEFINED HERE.
  ! It is allocated in minresqlptest* and used by Aprod* or Msolve.
  ! By default these variables are private.

  real(dp),    allocatable  :: d(:)     ! Defines diagonal matrix D.
  real(dp)                  :: Ashift   ! Shift diagonal elements of D in  Msolve.
  real(dp)                  :: Mpert    ! Perturbation to D in Msolve
                                        ! to avoid having an exact preconditioner.

  integer(ip)               :: nnz      ! These ones are used by the
  integer(ip), allocatable  :: indx(:)  ! Matrix Market Aprod routines
  integer(ip), allocatable  :: jndx(:)  ! AprodMtxCDS, AprodMtxCPS, AprodMtxCRS
  real(dp)   , allocatable  :: dval(:)  !
  integer(ip), allocatable  :: ival(:)  !
  real(sp)   , allocatable  :: rval(:)  !
  complex(dp), allocatable  :: cval(:)  !

contains

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine Aprod (n,x,y)

    integer(ip), intent(in)  :: n
    real(dp),    intent(in)  :: x(n)
    real(dp),    intent(out) :: y(n)

    !-------------------------------------------------------------------
    ! Aprod  computes y = A*x for some matrix A.
    ! This is a simple example for testing MINRESQLP.
    !-------------------------------------------------------------------

    integer(ip) :: i

    do i = 1, n
       y(i) = d(i)*x(i)
    end do

  end subroutine Aprod

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine Msolve(n,x,y)

    integer(ip), intent(in)  :: n
    real(dp),    intent(in)  :: x(n)
    real(dp),    intent(out) :: y(n)

    !-------------------------------------------------------------------
    ! Msolve solves M*y = x for some symmetric positive-definite matrix M.
    ! This is a simple example for testing MINRESQLP.
    ! Ashift will be the same as shift in MINRESQLP.
    !
    ! If Mpert = 0, the preconditioner will be exact, so
    ! MINRESQLP should require either one or two iterations,
    ! depending on whether (A - shift*I) is positive definite or not.
    !
    ! If Mpert is nonzero, somewhat more iterations will be required.
    !-------------------------------------------------------------------

    intrinsic   :: abs, mod
    integer(ip) :: i
    real(dp)    :: di

    do i = 1, n
       di   = abs( d(i) - Ashift )
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

  subroutine AprodMtxCDS (n,x,y)

    integer(ip), intent(in)  :: n
    real(dp),    intent(in)  :: x(n)
    real(dp),    intent(out) :: y(n)

    !-------------------------------------------------------------------
    ! AprodMtxCDS  computes y = A*x for some symmetric matrix A
    ! stored in Matrix Market CDS format (coordinate double symmetric).
    ! Only subdiagonal and diagonal elements are in (indx, jndx, dval).
    !-------------------------------------------------------------------

    integer(ip)  :: i, j, k
    real(dp)     :: d

    y(1:n) = zero

    do k = 1, nnz
       i = indx(k)
       j = jndx(k)
       d = dval(k)
!      d = rval(k)

       if (i > j) then          ! d = subdiagonal
          y(i) = y(i) + d*x(j)
          y(j) = y(j) + d*x(i)
       else                     ! i = j, d = diagonal
          y(i) = y(i) + d*x(i)
       end if

!      if (k <= 10) write(*,*) '  ', i, ' ', j, ' ', d,  ' ', y(i)
    end do

  end subroutine AprodMtxCDS

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine AprodMtxCPS (n,x,y)

    integer(ip), intent(in)  :: n
    real(dp),    intent(in)  :: x(n)
    real(dp),    intent(out) :: y(n)

    !-------------------------------------------------------------------
    ! AprodMtxCPS  computes y = A*x for some symmetric matrix A
    ! stored in Matrix Market CPS format (coordinate pattern symmetric).
    ! Only subdiagonal and diagonal elements are in (indx, jndx).
    ! The matrix entries are assumed to be 1.0.
    ! 12 Apr 2013: AprodMtxCPS derived from AprodMtxCRS.
    !-------------------------------------------------------------------

    integer(ip)  :: i, j, k

    y(1:n) = zero

    do k = 1, nnz
       i = indx(k)
       j = jndx(k)

       if (i > j) then          ! subdiagonal
          y(i) = y(i) + x(j)
          y(j) = y(j) + x(i)
       else                     ! i = j, diagonal
          y(i) = y(i) + x(i)
       end if

       ! if (k <= 10) write(*,*) '  ', i, ' ', j, ' ', y(i)
    end do

  end subroutine AprodMtxCPS

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine AprodMtxCRS (n,x,y)

    integer(ip), intent(in)  :: n
    real(dp),    intent(in)  :: x(n)
    real(dp),    intent(out) :: y(n)

    !-------------------------------------------------------------------
    ! AprodMtxCRS  computes y = A*x for some symmetric matrix A
    ! stored in Matrix Market CRS format (coordinate real symmetric).
    ! Only subdiagonal and diagonal elements are in (indx, jndx, rval).
    !-------------------------------------------------------------------

    integer(ip)  :: i, j, k
    real(dp)     :: d

    y(1:n) = zero

    do k = 1, nnz
       i = indx(k)
       j = jndx(k)
!      d = dval(k)
       d = rval(k)

       if (i > j) then          ! d = subdiagonal
          y(i) = y(i) + d*x(j)
          y(j) = y(j) + d*x(i)
       else                     ! i = j, d = diagonal
          y(i) = y(i) + d*x(i)
       end if

!      if (k <= 10) write(*,*) '  ', i, ' ', j, ' ', d,  ' ', y(i)
    end do

  end subroutine AprodMtxCRS

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine minresqlptest( n, precon, shift, pertM, sing, consis, default, nout )

    integer(ip), intent(in)    :: n, nout
    logical,     intent(in)    :: precon, sing, consis
    real(dp),    intent(in)    :: shift, pertM
    logical,     intent(in), optional :: default

    !-------------------------------------------------------------------
    ! minresqlptest solves sets up and solves a system (A - shift*I)x = b,
    ! using Aprod to define A and Msolve to define a preconditioner.
    !-------------------------------------------------------------------

    intrinsic :: real, present

    ! Local arrays and variables
    real(dp)    :: b(n), r1(n), w(n), x(n), xtrue(n), y(n)
    logical     :: disable, use_default
    integer(ip) :: j, itnlim, istop, itn, nprint
    real(dp)    :: Anorm, Acond, Arnorm, rnorm, rtol, r1norm, xnorm
    real(dp)    :: enorm, etol, wnorm, xnormtrue
    real(dp)    :: maxxnorm, TranCond, Acondlim

    character(len=*), parameter ::  headerStr =                                   &
       "(// '-----------------------------------------------------'"    //        &
       "  / 'Test of  MINRESQLP.'"                                      //        &
       "  / '-----------------------------------------------------'"    //        &
       "  / 'shift  =', f12.4, 6x, 'pertM  =', f12.4)"
    character(len=*), parameter ::  footerStr1 =                                  &
       "(// '  minresqlp appears to be successful.  n =', i7, '  Itns =', i7," // &
       "'  Relative error in x =', 1p, e8.1)"
    character(len=*), parameter ::  footerStr2 =                                  &
       "(// '  minresqlp appears to have failed.    n =', i7, '  Itns =', i7," // &
       "'  Relative error in x =', 1p, e8.1)"
    character(len=*), parameter ::  debugStr1 =                                   &
       "(// 'Final residual =', 1p, e8.1)"
    character(len=*), parameter ::  debugStr2 =                                   &
       "(// 'Solution  x', 1p, 4(i6, e14.6))"

    write(nout, headerStr) shift, pertM

    allocate( d(n) )           ! Array used in Aprod and Msolve to define A and M.
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
      xtrue(n-1:n) = zero
    end if

    call Aprod (n,xtrue,b)   ! b = A*xtrue
    b      = b - shift*xtrue ! Now b = (A - shift*I)*xtrue
    if (.not. consis) then
       b(n-1:n) = one        ! b not in the range of A
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

    disable  = .false.          ! Set other parameters and solve.
    itnlim   = n*3
    rtol     = 1.0e-12_dp
    maxxnorm = 1.0e+8_dp
    TranCond = 1.0e+7_dp
    Acondlim = 1.0e+15_dp

    if (prcsn == 6) then                     ! single precision is used
        rtol     = 1e-5_dp                   ! use a bigger rtol
        maxxnorm = 1.0e+4_dp                 ! use a smaller maxxnorm
        TranCond = 1.0e+3_dp
        Acondlim = 1.0e+7_dp
    end if

    if (debug) then
       write(*,*)
       write(*,*)  'n = ', n, ', precon = ', precon, ', shift = ', shift, ', pertM = ', pertM, &
                   ', sing = ', sing, ', consis = ', consis, ', nout = ', nout
       write(*,*)
       write(*,*)  'itnlim = ', itnlim, ', nout = ', nout,                &
                   ', maxxnorm = ', maxxnorm, ', TranCond = ', TranCond, ', Acondlim = ', Acondlim
    end if

    if (present(default)) then
        use_default = default
    else
        use_default = .false.
    end if

    if (use_default) then
       call MINRESQLP( n=n, Aprod=Aprod, b=b, shift=shift, x=x, nout=nout, itn=itn )
    else
       if (precon) then
         call MINRESQLP( n, Aprod, b, shift, Msolve, disable,                  &
                       nout, itnlim, rtol, maxxnorm, trancond, Acondlim,       &
                       x, istop, itn, rnorm, Arnorm, xnorm, Anorm, Acond )
       else
         call MINRESQLP( n=n, Aprod=Aprod, b=b, shift=shift, disable=disable,   &
                         nout=nout, itnlim=itnlim, rtol=rtol, maxxnorm=maxxnorm,&
                         trancond=trancond, Acondlim=Acondlim, x=x, istop=istop,&
                         itn=itn, rnorm=rnorm, Arnorm=Arnorm, xnorm=xnorm,      &
                         Anorm=Anorm, Acond=Acond )
       end if
    end if

    if (debug) then
       call Aprod (n,x,y)       ! y = A*x
       r1     = b - y + shift*x ! Final residual r1 = b - (A - shift*I)*x.
       r1norm = dnrm2(n,r1,1)
       write(nout,debugStr1) r1norm

       nprint = min(n,20)
       write(nout,*) ' '
       write(nout,*) 'Some of x'
       write(nout,debugStr2) (j, x(j), j=1,nprint)
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

    w      = x - xtrue            ! Print a clue about whether
    wnorm  = dnrm2(n,w,1)         ! the solution looks OK.
    xnormtrue = dnrm2(n,xtrue,1)
    enorm  = wnorm/xnormtrue

    if (use_default) then
       etol   = 1e-6
    else
       etol   = eps*n*Acond*10
    end if

    if (prcsn == 6) then                     ! single precision is used
        etol = 1e-3_dp                       ! use a bigger etol
    end if

    if (enorm <= etol) then
       write(nout, footerStr1) n, itn, enorm
    else
       write(nout, footerStr2) n, itn, enorm
    end if
    !write(nout, *)  etol, use_default
    deallocate(d)                           ! Free work array

  end subroutine minresqlptest

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine minresqlptestMtxCDS(input_file, consis, nout, tol)

    character(80),  intent(in)     :: input_file
    logical,        intent(in)     :: consis
    integer(ip),    intent(in)     :: nout
    real(dp),       intent(in)     :: tol

    !-------------------------------------------------------------------
    ! 29 Oct 2012: Use xcheck to check computed x from MINRESQLP.
    ! 02 Jan 2013: Print n in the " minresqlp  appears to be ..." message
    !              to help identify the problem.
    ! 21 Apr 2013: minresqlptestMtxCDS created from minresqlpmtxtest.
    !-------------------------------------------------------------------

    intrinsic      :: real

    integer(ip)    :: input_unit, nrow, ncol
    character(14)  :: id
    character(10)  :: rep
    character( 6)  :: type
    character(7)   :: field
    character(19)  :: symm

    real(dp), allocatable  :: b(:), x(:), r1(:), w(:)

    logical        :: disable
    integer(ip)    :: n, j, itn, itnlim, istop, inform
    real(dp)       :: shift, Anorm, Acond, Arnorm, rnorm, rtol, xnorm
    real(dp)       :: maxxnorm, TranCond, Acondlim
    real(dp)       :: test1, test2, tolcheck

    character(len=*), parameter :: headerStr =             &
       "(// '---------------------------------------'"  // &
       "  / ' Test of MINRESQLP on an MM CDS matrix '"  // &
       "  / '---------------------------------------')"
    character(len=*), parameter :: footerStr1 = &
       "(// '  minresqlp appears to be successful.  n =', i7, '  Itns =', i7," // &
       "'  test(r) =', 1p, e9.2, '  test(Ar) =', 1p, e9.2)"
    character(len=*), parameter :: footerStr2 = &
       "(// '  minresqlp appears to have failed.    n =', i7, '  Itns =', i7," // &
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
    allocate( rval(nnz) )   ! CDS
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

    do j = 1, n
       b(j) = one    ! real(j,dp)/real(n,dp)
       x(j) = real(j,dp)/n
    end do

    !write(nout,*) (j, b(j), j=1,5)  ! Print some of the solution

    !call AprodMtxCDS (n,b,x)   ! b = A*xtrue
    !write(nout,debugStr) (j, x(j), j=1,10)  ! Print some of the solution
    !call AprodMtxCDS (n,x,b)
    !write(nout,debugStr) (j, b(j), j=1,10)  ! Print some of the solution

    if (debug) then
      write(nout,*) 'consis     = ', consis
      write(nout,*) 'input_file = ', trim(input_file)
      write(nout,*) 'n = ', n, '  nnz = ', nnz
    end if

    disable  = .false.          ! Set other parameters and solve.
    itnlim   = n*3
    rtol     = tol
    maxxnorm = 1.0e+8_dp
    TranCond = 1.0e+8_dp
    Acondlim = 1.0e+15_dp
    shift    = zero

    if (prcsn == 6) then                     ! single precision is used
        rtol     = 1e-5_dp                   ! use a bigger rtol
        maxxnorm = 1.0e+4_dp                 ! use a smaller maxxnorm
        TranCond = 1.0e+3_dp
        Acondlim = 1.0e+7_dp
    end if

    if (consis) then
       call AprodMtxCDS(n,x,b)   ! b = A*x
       if (debug) then
         write(nout,*) ' '
         write(nout,*) 'norm(b) =', dnrm2(n,b,1)
         write(nout,*) 'Some of the x defining b'
         do j = 1, min(n,5)
            write(nout,*) j, x(j)
         end do
       end if
    end if

    if (debug) then
      write(*,*) 'Some of b'
      do j = 1, min(n,5)
         write(nout,*) j, b(j)
      end do
    end if

    if (debug) then
       write(*,*)
       write(*,*)  'n = ', n, ', shift = ', shift,  ', consis = ', consis, ', nout = ', nout
       write(*,*)
       write(*,*)  'itnlim = ', itnlim, ', nout = ', nout,                &
                   ', maxxnorm = ', maxxnorm, ', TranCond = ', TranCond, ', Acondlim = ', Acondlim
    end if


    call MINRESQLP( n=n, Aprod=AprodMtxCDS, b=b, shift=shift, disable=disable,&
                    nout=nout, itnlim=itnlim, rtol=rtol, maxxnorm=maxxnorm,    &
                    trancond=trancond, Acondlim=Acondlim, x=x, istop=istop,   &
                    itn=itn, rnorm=rnorm, Arnorm=Arnorm, xnorm=xnorm,         &
                    Anorm=Anorm, Acond=Acond )

    tolcheck = 10_dp*tol
    call xcheck( n, AprodMtxCDS, b, shift, x, Anorm, tolcheck, nout, &
                 test1, test2, inform )

    if (inform <= 2) then
       write(nout, footerStr1) n, itn, test1, test2
    else
       write(nout, footerStr2) n, itn, test1, test2
    end if

    if (debug) then
      write(nout,*) ' '
      write(nout,*) 'norm(x) =', dnrm2(n,x,1)
      write(nout,*) 'Some of the computed x'
      do j = 1, min(n,5)
         write(nout,*) x(j)
      end do
    end if

    deallocate( indx, jndx, ival, rval, dval, cval )
    deallocate( b, x, r1, w )

  end subroutine minresqlptestMtxCDS

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine minresqlptestMtxCPS(input_file, consis, nout, tol)

    character(80),  intent(in)     :: input_file
    logical,        intent(in)     :: consis
    integer(ip),    intent(in)     :: nout
    real(dp),       intent(in)     :: tol

    !-------------------------------------------------------------------
    ! 29 Oct 2012: Use xcheck to check computed x from MINRESQLP.
    ! 02 Jan 2013: Print n in the " minresqlp  appears to be ..." message
    !              to help identify the problem.
    ! 21 Apr 2013: minresqlptestMtxCPS created from minresqlpmtxtest.
    !-------------------------------------------------------------------

    intrinsic      :: real

    integer(ip)    :: input_unit, nrow, ncol
    character(14)  :: id
    character(10)  :: rep
    character( 6)  :: type
    character(7)   :: field
    character(19)  :: symm

    real(dp), allocatable  :: b(:), x(:), r1(:), w(:)

    logical        :: disable
    integer(ip)    :: n, j, itn, itnlim, istop, inform
    real(dp)       :: shift, Anorm, Acond, Arnorm, rnorm, rtol, xnorm
    real(dp)       :: maxxnorm, TranCond, Acondlim
    real(dp)       :: test1, test2, tolcheck

    character(len=*), parameter :: headerStr =             &
       "(// '---------------------------------------'"  // &
       "  / ' Test of MINRESQLP on an MM CPS matrix '"  // &
       "  / '---------------------------------------')"
    character(len=*), parameter :: footerStr1 = &
       "(// '  minresqlp appears to be successful.  n =', i7, '  Itns =', i7," // &
       "'  test(r) =', 1p, e9.2, '  test(Ar) =', 1p, e9.2)"
    character(len=*), parameter :: footerStr2 = &
       "(// '  minresqlp appears to have failed.    n =', i7, '  Itns =', i7," // &
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
    allocate( rval(nnz) )   ! CPS
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

    do j = 1, n
       b(j) = one    ! real(j,dp)/real(n,dp)
       x(j) = real(j,dp)/n
    end do

    !call AprodMtxCPS (n,b,x)   ! b = A*xtrue
    !write(nout,debugStr) (j, x(j), j=1,10)  ! Print some of the solution
    !call AprodMtxCPS (n,x,b)
    !write(nout,debugStr) (j, b(j), j=1,10)  ! Print some of the rhs

    if (debug) then
      write(nout,*) 'consis     = ', consis
      write(nout,*) 'input_file = ', trim(input_file)
      write(nout,*) 'n = ', n, '  nnz = ', nnz
    end if

    disable  = .false.          ! Set other parameters and solve.
    itnlim   = n*3
    rtol     = tol
    maxxnorm = 1.0e+8_dp
    TranCond = 1.0e+8_dp
    Acondlim = 1.0e+15_dp
    shift    = zero

    if (prcsn == 6) then                     ! single precision is used
        rtol     = 1e-5_dp                   ! use a bigger rtol
        maxxnorm = 1.0e+4_dp                 ! use a smaller maxxnorm
        TranCond = 1.0e+3_dp
        Acondlim = 1.0e+7_dp
    end if

    if (consis) then
       call AprodMtxCPS (n,x,b)   ! b = A*x
       if (debug) then
         write(nout,*) ' '
         write(nout,*) 'norm(b) =', dnrm2(n,b,1)
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
       write(*,*)  'itnlim = ', itnlim, ', nout = ', nout,                &
                   ', maxxnorm = ', maxxnorm, ', TranCond = ', TranCond, ', Acondlim = ', Acondlim
    end if

 !   call MINRESQLP( n, AprodMtxCPS, b, shift, Msolve, disable,           &
 !                   nout, itnlim, rtol, maxxnorm, trancond, Acondlim,    &
 !                   x, istop, itn, rnorm, Arnorm, xnorm, Anorm, Acond )

    call MINRESQLP( n=n, Aprod=AprodMtxCPS, b=b, shift=shift, disable=disable,&
                    nout=nout, itnlim=itnlim, rtol=rtol, maxxnorm=maxxnorm,   &
                    trancond=trancond, Acondlim=Acondlim, x=x, istop=istop,   &
                    itn=itn, rnorm=rnorm, Arnorm=Arnorm, xnorm=xnorm,         &
                    Anorm=Anorm, Acond=Acond )

    tolcheck = 10_dp*tol
    call xcheck( n, AprodMtxCPS, b, shift, x, Anorm, tolcheck, nout, &
                 test1, test2, inform )

    if (inform <= 2) then
       write(nout, footerStr1) n, itn, test1, test2
    else
       write(nout, footerStr2) n, itn, test1, test2
    end if

    if (debug) then
      write(nout,*) ' '
      write(nout,*) 'norm(x) =', dnrm2(n,x,1)
      write(nout,*) 'Some of the computed x'
      do j = 1, min(n,5)
         write(nout,*) j, x(j)
      end do
    end if

    deallocate( indx, jndx, ival, rval, dval, cval )
    deallocate( b, x, r1, w )

  end subroutine minresqlptestMtxCPS

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine minresqlptestMtxCRS(input_file, consis, nout, tol)

    character(80),  intent(in)     :: input_file
    logical,        intent(in)     :: consis
    integer(ip),    intent(in)     :: nout
    real(dp),       intent(in)     :: tol

    !-------------------------------------------------------------------
    ! 29 Oct 2012: Use xcheck to check computed x from MINRESQLP.
    ! 02 Jan 2013: Print n in the " minresqlp  appears to be ..." message
    !              to help identify the problem.
    !-------------------------------------------------------------------

    intrinsic      :: real

    integer(ip)    :: input_unit, nrow, ncol
    character(14)  :: id
    character(10)  :: rep
    character( 6)  :: type
    character(7)   :: field
    character(19)  :: symm

    real(dp), allocatable  :: b(:), x(:), r1(:), w(:)

    logical        :: disable
    integer(ip)    :: n, j, itn, itnlim, istop, inform
    real(dp)       :: shift, Anorm, Acond, Arnorm, rnorm, rtol, xnorm
    real(dp)       :: maxxnorm, TranCond, Acondlim
    real(dp)       :: test1, test2, tolcheck

    character(len=*), parameter :: headerStr =             &
       "(// '---------------------------------------'"  // &
       "  / ' Test of MINRESQLP on an MM CRS matrix '"  // &
       "  / '---------------------------------------')"
    character(len=*), parameter :: footerStr1 = &
       "(// '  minresqlp appears to be successful.  n =', i7, '  Itns =', i7," // &
       "'  test(r) =', 1p, e9.2, '  test(Ar) =', 1p, e9.2)"
    character(len=*), parameter :: footerStr2 = &
       "(// '  minresqlp appears to have failed.    n =', i7, '  Itns =', i7," // &
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

    do j = 1, n
       b(j) = one    ! real(j,dp)/real(n,dp)
       x(j) = real(j,dp)/n
    end do

    !write(nout,*) (j, b(j), j=1,5)  ! Print some of the solution

    !call AprodMtxCRS (n,b,x)   ! b = A*xtrue
    !write(nout,debugStr) (j, x(j), j=1,10)  ! Print some of the solution
    !call AprodMtxCRS (n,x,b)
    !write(nout,debugStr) (j, b(j), j=1,10)  ! Print some of the solution

    if (debug) then
      write(nout,*) 'consis     = ', consis
      write(nout,*) 'input_file = ', trim(input_file)
      write(nout,*) 'n = ', n, '  nnz = ', nnz
    end if

    disable  = .false.          ! Set other parameters and solve.
    itnlim   = n*3
    rtol     = tol
    maxxnorm = 1.0e+8_dp
    TranCond = 1.0e+8_dp
    Acondlim = 1.0e+15_dp
    shift    = zero

    if (prcsn == 6) then                     ! single precision is used
        rtol     = 1e-5_dp                   ! use a bigger rtol
        maxxnorm = 1.0e+4_dp                 ! use a smaller maxxnorm
        TranCond = 1.0e+3_dp
        Acondlim = 1.0e+7_dp
    end if

    if (consis) then
       call AprodMtxCRS (n,x,b)   ! b = A*x
       if (debug) then
         write(nout,*) ' '
         write(nout,*) 'norm(b) =', dnrm2(n,b,1)
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
       write(*,*)  'itnlim = ', itnlim, ', nout = ', nout,                &
                   ', maxxnorm = ', maxxnorm, ', TranCond = ', TranCond, ', Acondlim = ', Acondlim
    end if


    call MINRESQLP( n=n, Aprod=AprodMtxCRS, b=b, shift=shift, disable=disable,&
                    nout=nout, itnlim=itnlim, rtol=rtol, maxxnorm=maxxnorm,   &
                    trancond=trancond, Acondlim=Acondlim, x=x, istop=istop,   &
                    itn=itn, rnorm=rnorm, Arnorm=Arnorm, xnorm=xnorm,         &
                    Anorm=Anorm, Acond=Acond )

    tolcheck = 10_dp*tol
    call xcheck( n, AprodMtxCRS, b, shift, x, Anorm, tolcheck, nout, &
                 test1, test2, inform )

    if (inform <= 2) then
       write(nout, footerStr1) n, itn, test1, test2
    else
       write(nout, footerStr2) n, itn, test1, test2
    end if

    if (debug) then
      write(nout,*) ' '
      write(nout,*) 'norm(x) =', dnrm2(n,x,1)
      write(nout,*) 'Some of the computed x'
      do j = 1, min(n,5)
         write(nout,*) j, x(j)
      end do
    end if

    deallocate( indx, jndx, ival, rval, dval, cval )
    deallocate( b, x, r1, w )

  end subroutine minresqlptestMtxCRS

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine symorthotest( a, b, c_true, s_true, r_true, nout, tol )
    real(dp),    intent(in)           :: a, b, c_true, s_true, r_true
    real(dp),    intent(in), optional :: tol
    integer(ip), intent(in)           :: nout

    !-------------------------------------------------------------------
    ! 20 Aug 2012: First version of symorthotest.
    !-------------------------------------------------------------------

    real(dp)                    :: c, s, r, norm_diff
    real(dp)                    :: out(3), expected_out(3), relTol
    character(len=*), parameter :: footerStr1 = &
       "(// '  symortho  appears to be successful.  Relative error in [c,s,r] =', 1p, e8.1)"
    character(len=*), parameter :: footerStr2 = &
       "(// '  symortho  appears to have failed.    Relative error in [c,s,r] =', 1p, e8.1)"

    call symortho( a, b, c, s, r )

    out          = (/    c  , s     , r      /)
    expected_out = (/ c_true, s_true, r_true /)
    norm_diff    = dnrm2( 3, out - expected_out, 1 ) / dnrm2( 3, expected_out, 1 )

    if (.not. present(tol)) then
       relTol = eps
    else
       relTol = tol
    end if

    write(nout,*)  ' '
    write(nout,*)  '-----------------------------------------------------'
    write(nout,*)  'Test of  SYMORTHO.'
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

  end subroutine symorthotest

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! subroutine   xcheck
!
! xcheck tests if a given x seems to be a solution of symmetric system
! Ax = b or Ax ~= b.
!
! 29 Oct 2012: derived from lsmrCheckModule.f90.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine xcheck( n, Aprod, b, shift, x, Anorm, tol, nout, &
                     test1, test2, inform )

    integer(ip), intent(in)    :: n        ! No. of rows and cols of A
    integer(ip), intent(in)    :: nout     ! Output file number
    integer(ip), intent(out)   :: inform   ! = 0 if b = 0 and x = 0.
                                           ! = 1 or 2 if x seems to
                                           !   solve systems 1 or 2 below;
                                           ! = 4 otherwise.
    real(dp),    intent(in)    :: Anorm    ! An estimate of norm(A - shift*I)
                                           ! provided by MINRESQLP.
    real(dp),    intent(in)    :: tol      ! tolerance for judging residuals.
                                           ! Typically the tol used for computing x.
    real(dp),    intent(in)    :: shift    ! Often zero.
    real(dp),    intent(in)    :: b(n)     ! The right-hand side of (A - shift*I)x ~= b.
    real(dp),    intent(in)    :: x(n)     ! The given solution estimate.
    real(dp),    intent(out)   :: test1    ! Should be small for consistent systems.
    real(dp),    intent(out)   :: test2    ! Should be small for LS systems.

    interface
       subroutine Aprod(n,x,y)             ! y := A*x
         use minresqlpDataModule, only : dp, ip
         integer(ip), intent(in)    :: n
         real(dp),    intent(in)    :: x(n)
         real(dp),    intent(out)   :: y(n)
       end subroutine Aprod
    end interface

    !-------------------------------------------------------------------
    ! One-liner: xcheck tests if x solves symmetric (A-shift*I)x=b.
    !
    ! Purpose:   xcheck computes residuals and norms associated with
    ! the vector x and the problem solved by MINRESQLP.
    ! It determines whether x seems to be a solution to either of the
    ! systems:  1. (A - shift*I)x = b
    !           2. min ||(A - shift*I)x - b||
    !
    ! History:
    ! 29 Oct 2012: xcheck derived from xcheck in lsmrCheckModule.f90.
    ! 02 Jan 2013: Print Anorm, xnorm, bnorm.
    !-------------------------------------------------------------------

    intrinsic           :: epsilon, max

    ! Local variables and arrays
    real(dp)            :: r(n), v(n)
    real(dp)            :: bnorm, eps, rnorm, sigma, tol2, xnorm

    eps    = epsilon(eps)
    tol2   = max( tol, eps )

    call Aprod(n,x,r)        ! r = Ax
    r      = b - r + shift*x ! r = b - (A - shift*I)x = b - Ax + shift*x
    call Aprod(n,r,v)        ! v = Ar
    v      = v - shift*r     ! v = (A - shift*I)r = Ar - shift*r

    bnorm  = dnrm2 (n,b,1)   ! Compute the norms of b, x, r, v.
    xnorm  = dnrm2 (n,x,1)
    rnorm  = dnrm2 (n,r,1)
    sigma  = dnrm2 (n,v,1)

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
      // ' Enter xcheck.  Does x solve Ax = b?  where A is really (A - shift*I)' &
      /  '    shift     =', e10.3  &
      /  '    norm(A)   =', e10.3  &
      /  '    norm(x)   =', e10.3  &
      /  '    norm(b)   =', e10.3  &
      /  '    norm(r)   =', e15.8  &
      /  '    norm(Ar)  =', e10.3)
 3000 format(1p                    &
      /  '    inform    =', i2     &
      /  '    tol       =', e10.3  &
      /  '    test1     =', e10.3, ' (Ax = b)           rnorm/(Anorm*xnorm + bnorm)' &
      /  '    test2     =', e10.3, ' (least-squares)    Arnorm/(Anorm*rnorm)')

  end subroutine xcheck

end module minresqlpTestModule
