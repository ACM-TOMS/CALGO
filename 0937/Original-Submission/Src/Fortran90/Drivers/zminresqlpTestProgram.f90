!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! File zminresqlpTestProgram.f90
!
!    zminresqlpTestProgram
!
!
! This program calls zminresqlptest(...) and optionally zminresqlptestMtxCCH(...) in
! zminresqlpTestModule to generate a series of test problems
! Ax = b or Ax ~= b and solve them with ZMINRESQLP.
! The matrix A is n x n.  It is defined by routines in zminresqlpTestModule.
!
! There are also optional tests for complex SYMORTHO.
!
! To turn on/off optional tests for SYMORTHO and MINRESQLP using files
! in matrix-market format, change the values of local logical constants
! testSymortho and testMtx in this file.
!
! To compile this program:
!    make -f zMakefile
! To run this program:
!    ./zminresqlptest
! Outputs from up to 63 unit tests is written to zMINRESQLP.txt.
! To check results, search for the word "appear" by running:
!    grep appear zMINRESQLP.txt | cat -n
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
! 12 Jul 2011: Created complex version zminresqlpTestProgram.f90
!              from real version minresqlpTestProgram.f90.
! 20 Aug 2012: Added calls to two subroutines minresqlpmtxtest and symorthotest
!              for testing SYMORTHO and MINRESQLP on singular matrices
!              in Matrix Market Format.
! 15 Sep 2012: Matrix dwg961a.mtx from SJSU collection of singular matrices
!              changed from CCS format (coordinate complex symmtric)
!                        to CCH format (coordinate complex hermitian)
!              to provide example DataMtx/CCH/dwg961aCCH.mtx.
! 15 Sep 2012: Matrix qc324.mtx from UFL collection of matrices
!              changed from CCS format (coordinate complex symmtric)
!                        to CCH format (coordinate complex hermitian)
!              to provide example DataMtx/CCH/qc324CCH.mtx.
! 30 Oct 2012: shift is now real(d), not complex(dp).
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

program zminresqlpTestProgram

  use   zminresqlpDataModule, only : zzero, zone, dp, ip, prcsn, zero, eps, one, &
                                     realmin, realmax, testSymortho, testMtx
  use   zminresqlpTestModule, only : zminresqlptest,       &
                                     zminresqlptest2,      &
                                     zminresqlptestMtxCCH, &
                                     zsymorthotest

  implicit none

  intrinsic :: trim, sqrt, conjg, abs

  ! Local variables
  logical                :: normal, precon, sing, consis, default
  integer(ip)            :: n, nout, ios, i, j
  complex(dp)            :: a, b, s_true, r_true, r_val, BIG, SMALL
  real(dp)               :: c_true, shift, pertM, tol
  integer(ip), parameter :: nCCH = 3
  character(len=*), parameter :: pathCCH = './DataMtx/CCH/'
  character(len=*), parameter :: input_files(nCCH) =  &
      (/               &
         'toy5      '  &    ! n =    5
       , 'qc324CCH  '  &    ! n =  324
       , 'dwg961aCCH'  &    ! n =  961
       /)
  character(80)          :: input_file, output_file

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  write(*,*) "integer and real precision:"
  write(*,*) "ip =", ip
  write(*,*) "dp =", dp

  nout   = 10
  output_file = 'zMINRESQLP.txt'
  open(nout, file=output_file, status='unknown', iostat=ios)

  if (ios /= 0) then
    write(*,*) "Error opening file ", trim(output_file)
  end if

  !------------------------------------------------------------------
  ! tests for zsymortho
  !------------------------------------------------------------------
  if (testSymortho) then
    a      = zzero
    b      = zzero
    c_true = one
    s_true = zzero
    r_true = zzero
    call zsymorthotest ( a, b, c_true, s_true, r_true, nout )            ! Ex 01

    a      = zzero
    b      = 2_dp*zone
    c_true = zero
    s_true = zone
    r_true = 2_dp*zone
    call zsymorthotest ( a, b, c_true, s_true, r_true, nout )            ! Ex 02

    a      = zone *     6_dp
    b      = zone * (  -8_dp)
    c_true =  one *   0.6_dp
    s_true = zone * (-0.8_dp)
    r_true = zone * 10_dp
    call zsymorthotest ( a, b, c_true, s_true, r_true, nout )            ! Ex 03

    a      = zone *   8_dp
    b      = zone *   6_dp
    c_true =  one * 0.8_dp
    s_true = zone * 0.6_dp
    r_true = zone *  10_dp
    call zsymorthotest ( a, b, c_true, s_true, r_true, nout )            ! Ex 04

    BIG    = zone * realmax / 1000.0_dp
    a      =  BIG * (-3_dp)
    b      =  BIG * (-4_dp)
    r_true =  BIG * (-5_dp)
    c_true = abs(a) / abs(r_true)
    s_true = b / r_true
    call zsymorthotest ( a, b, c_true, s_true, r_true, nout )            ! Ex 05

    SMALL  = zone * realmin * 1000.0_dp
    a      = SMALL * (-3_dp)
    b      = SMALL * (-4_dp)
    r_true = SMALL * (-5_dp)
    c_true = abs(a) / abs(r_true)
    s_true = b / r_true
    call zsymorthotest ( a, b, c_true, s_true, r_true, nout )            ! Ex 06

    a      =  BIG * (-3_dp)
    b      = zone * (-4_dp)
    r_true =  BIG * (-3_dp)
    c_true = abs(a) / abs(r_true)
    s_true = b / r_true
    call zsymorthotest ( a, b, c_true, s_true, r_true, nout )            ! Ex 07

    a      = zone * (-4_dp)
    b      =  BIG * (-3_dp)
    r_true =  BIG * (-3_dp)
    c_true = abs(a) / abs(r_true)
    s_true = b / r_true
    call zsymorthotest ( a, b, c_true, s_true, r_true, nout )            ! Ex 08

    a      = zone  * (-3_dp)
    b      = SMALL * (-4_dp)
    r_true = zone  * (-3_dp)
    c_true = abs(a) / abs(r_true)
    s_true = b / r_true
    call zsymorthotest ( a, b, c_true, s_true, r_true, nout )            ! Ex 09

    a      = SMALL * (-3_dp)
    b      = zone  * (-4_dp)
    r_true = zone  * (-4_dp)
    c_true = abs(a) / abs(r_true)
    s_true = b / r_true
    call zsymorthotest ( a, b, c_true, s_true, r_true, nout )            ! Ex 10

    ! complex numbers
    a      = cmplx(0.0_dp,1.0_dp,dp)
    b      = zzero
    r_true = a
    c_true = abs(a) / abs(r_true)
    s_true = b / r_true
    call zsymorthotest ( a, b, c_true, s_true, r_true, nout )            ! Ex 11

    a      = zzero
    b      = cmplx(0.0_dp,1.0_dp,dp)
    r_true = b
    c_true = abs(a) / abs(r_true)
    s_true = b / r_true
    call zsymorthotest ( a, b, c_true, s_true, r_true, nout )            ! Ex 12

    a      = cmplx(3.0_dp,4.0_dp,dp)
    b      = cmplx(4.0_dp,6.0_dp,dp)
    r_true = (a/abs(a)) * sqrt(abs(a)**2 + abs(b)**2)
    c_true = abs(a) / abs(r_true)
    s_true = conjg(b / r_true)
    call zsymorthotest ( a, b, c_true, s_true, r_true, nout )            ! Ex 13

    a      = cmplx(-3.0_dp,-4.0_dp,dp)
    b      = cmplx(-4.0_dp,-6.0_dp,dp)
    r_true = (a/abs(a)) * sqrt(abs(a)**2 + abs(b)**2)
    r_val  = r_true
    c_true = abs(a) / abs(r_true)
    s_true = conjg(b / r_true)
    call zsymorthotest ( a, b, c_true, s_true, r_true, nout )            ! Ex 14

    a      = BIG * cmplx(-3.0_dp,-4.0_dp,dp)
    b      = BIG * cmplx(-4.0_dp,-6.0_dp,dp)
    r_true = BIG * r_val
    c_true = abs(a) / abs(r_true)
    s_true = conjg(b / r_true)
    tol = 10 * eps
    call zsymorthotest ( a, b, c_true, s_true, r_true, nout, tol)        ! Ex 15

    a      = BIG * cmplx(-3.0_dp,-4.0_dp,dp)
    b      =       cmplx(-4.0_dp,-6.0_dp,dp)
    r_true = a
    c_true = abs(a) / abs(r_true)
    s_true = conjg(b / r_true)
    call zsymorthotest ( a, b, c_true, s_true, r_true, nout )            ! Ex 16

    a      =       cmplx(-4.0_dp,-6.0_dp,dp)
    b      = BIG * cmplx(-4.0_dp,-6.0_dp,dp)
    r_true = b
    c_true = abs(a) / abs(r_true)
    s_true = conjg(b / r_true)
    call zsymorthotest ( a, b, c_true, s_true, r_true, nout, tol )        ! Ex 17

    a      =         cmplx(-3.0_dp,-4.0_dp,dp)
    b      = SMALL * cmplx(-4.0_dp,-6.0_dp,dp)
    r_true = a
    c_true = abs(a) / abs(r_true)
    s_true = conjg(b / r_true)
    call zsymorthotest ( a, b, c_true, s_true, r_true, nout )            ! Ex 18

    a      = SMALL * cmplx(-3.0_dp,-4.0_dp,dp)
    b      =         cmplx(-4.0_dp,-6.0_dp,dp)
    r_true = (a/abs(a)) * sqrt(abs(a)**2 + abs(b)**2)
    c_true = abs(a) / abs(r_true)
    s_true = conjg(b / r_true)
    call zsymorthotest ( a, b, c_true, s_true, r_true, nout, tol)        ! Ex 19

    a      = SMALL * cmplx(-3.0_dp,-4.0_dp,dp)
    b      = SMALL * cmplx(-4.0_dp,-6.0_dp,dp)
    r_true = SMALL * r_val
    c_true = abs(a) / abs(r_true)
    s_true = conjg(b / r_true)
    call zsymorthotest ( a, b, c_true, s_true, r_true, nout, tol)        ! Ex 20


  end if

  !------------------------------------------------------------------
  ! Tests for ZMINRESQLP
  !------------------------------------------------------------------
  do i = 1, 2
     normal = .false. ! no preconditioning
     precon = .true.
     shift  = 0.01_dp
     pertM  = 0.1_dp
     sing   = .false.
     consis = .true.
     default = .false.
     if (i > 1) default = .true.

     write(nout,*) ' '
     write(nout,*) 'ZMINRESQLP tests with default = ', default

     ! Test the unlikely tiny cases that often trip us up.

     n      = 1
     call zminresqlptest( n, normal, zero,  zero, sing, consis, default, nout ) ! Ex 01
     call zminresqlptest( n, normal, shift, zero, sing, consis, default, nout ) ! Ex 02

     n      = 2
     call zminresqlptest( n, normal, zero,  zero, sing, consis, default, nout ) ! Ex 03
     call zminresqlptest( n, normal, shift, zero, sing, consis, default, nout ) ! Ex 04

     n      = 3
     call zminresqlptest( n, normal, zero,  zero, sing, consis, default, nout ) ! Ex 05
     call zminresqlptest( n, normal, shift, zero, sing, consis, default, nout ) ! Ex 06

     n      = 4
     call zminresqlptest( n, normal, zero,  zero, sing, consis, default, nout ) ! Ex 07
     call zminresqlptest( n, normal, shift, zero, sing, consis, default, nout ) ! Ex 08

     ! Test small positive-definite and indefinite systems
     ! without preconditioners.  ZMINRESQLP should take n iterations.

     n      = 50
     call zminresqlptest( n, normal, zero,  zero, sing, consis, default, nout ) ! Ex 09
     call zminresqlptest( n, normal, shift, zero, sing, consis, default, nout ) ! Ex 10

     ! Test small positive-definite and indefinite systems with
     ! exact preconditioners.  ZMINRESQLP should take about n iterations.

     n      = 2
     call zminresqlptest( n, precon, zero,  zero, sing, consis, default, nout ) ! Ex 11
     call zminresqlptest( n, precon, shift, zero, sing, consis, default, nout ) ! Ex 12

     n      = 50
     call zminresqlptest( n, precon, zero,  zero, sing, consis, default, nout ) ! Ex 13
     call zminresqlptest( n, precon, shift, zero, sing, consis, default, nout ) ! Ex 14

     ! pertM makes the preconditioners incorrect in n/10 entries.
     ! ZMINRESQLP should take about n/10 iterations.

     call zminresqlptest( n, precon, zero,  pertM, sing, consis, default, nout ) ! Ex 15
     call zminresqlptest( n, precon, shift, pertM, sing, consis, default, nout ) ! Ex 16

     ! Singular consistent test case.
     sing   = .true.
     consis = .true.

     n      = 4
     call zminresqlptest( n, normal, zero,  zero, sing, consis, default, nout )  ! Ex 17

     ! if shift = 0.25 or 0.5 and n =4, then the problem is inconsistent
     shift  = 0.3_dp
     call zminresqlptest( n, normal, shift, zero, sing, consis, default, nout )  ! Ex 18

     n      = 50
     call zminresqlptest( n, normal, zero,  zero, sing, consis, default, nout )  ! Ex 19

     ! Singular inconsistent test case.

     consis = .false.

     n      = 4
     call zminresqlptest( n, normal, zero, zero, sing, consis, default, nout ) ! Ex 20

     n      = 50
     call zminresqlptest( n, normal, zero, zero, sing, consis, default, nout ) ! Ex 21
  end do

  ! (singular) Hermitian inconsistent test case.
  n      = 5
  shift  = zero
  call zminresqlptest2( n, normal, shift, zero, nout )                       ! Ex 22

  !------------------------------------------------------------------
  ! Test ZMINRESQLP on Matrix Market files
  !------------------------------------------------------------------
  tol = 1e-7_dp       ! 8e-4_dp
  if (prcsn == 6) then                     ! single precision is used
    tol = 2.2e-3_dp   ! 3.7e-3_dp
  end if

  if (testMtx) then
     write(nout,*) ' '
     write(nout,*) 'ZMINRESQLP tests on MM CCH examples'

     do j = 1, 2
        consis = .true.
        if (j > 1) consis = .false.
        do i = 1, nCCH
           input_file = pathCCH//trim(input_files(i))//'.mtx'
           call zminresqlptestMtxCCH(input_file, consis, nout, tol)
        end do
     end do
  end if

  write(*,*)
  write(*,*) "Results are in output file  ", trim(output_file)
  write(*,*) "Search the file for 'appear'"
  write(*,*) "For example:    grep appear ", trim(output_file)
  write(*,*)

  close( unit = nout )
end program zminresqlpTestProgram
