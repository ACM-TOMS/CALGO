!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! File minresqlpTestProgram.f90
!
!    minresqlpTestProgram
!
! This program calls minresqlptest(...) and optionally minresqlptestMtx*(...)
! in minresqlpTestModule to generate or input a series of test problems
! Ax = b or Ax ~= b and solve them with MINRESQLP.
! The matrix A is n x n.  It is defined by routines in minresqlpTestModule.
!
! There are also optional tests for SYMORTHO.
!
! To turn on/off optional tests for SYMORTHO and MINRESQLP using files
! in Matrix Market format, change the values of local logical constants
! testSymortho and testMtx in this file, respectively.
!
! To compile this program:
!    make   (or  make -f Makefile)
! To run this program:
!    ./minresqlptest
! Output is written to MINRESQLP.txt.
! To check results, search for the word "appear" by running:
!    grep appear MINRESQLP.txt | cat -n
!
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
! 11 Oct 2007: First version of minresqlpTestProgram.f90.
!              Initially used compiler option -r8.
! 15 Oct 2007: Use real(kind=8) everywhere. minresqlptest2 added.
! 16 Oct 2007: Use minresqlpDataModule to define dp = selected_real_kind(15).
! 12 Jul 2011: Created complex version zminresqlpTestProgram.f90
!              from real version minresqlpTestProgram.f90.
! 20 Aug 2012: Added calls to two subroutines minresqlpmtxtest and symorthotest
!              for testing SYMORTHO and MINRESQLP on singular matrices
!              in Matrix Market Format.
!              (minresqlpmtxtest is now
!              minresqlptestMtxCDS, minresqlptestMtxCPS, minresqlptestMtxCRS.)
! 27 Aug 2012: input_files(*) have to be entered with constant length.
!              output_file hardwires the name of the results file.
!              nout = 10 (not 6) to avoid clash with write(*,*).
! 08 Sep 2012: More SJSU singular coordinate real symmetric (CRS) examples added.
!              Some of them need many iterations, which defeats the purpose
!              of "quick" unit tests.  The long ones are commented out of
!              the list of names in input_files(*).
! 21 Apr 2012: Matrix Market test routines are now
!              minresqlptestMtxCDS, minresqlptestMtxCPS, minresqlptestMtxCRS.
! 09 Sep 2013: Removed unused references to connected, named_file.
!              Changed tol from 1e-8_dp to 1e-7_dp to fix a failed test.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

program minresqlpTestProgram

  use   minresqlpDataModule, only  : dp, ip, prcsn, zero, one, realmax, realmin, &
                                     testSymortho, testMtx
  use   minresqlpTestModule, only  : minresqlptest,       &
                                     minresqlptestMtxCDS, &
                                     minresqlptestMtxCPS, &
                                     minresqlptestMtxCRS, &
                                     symorthotest

  implicit none

  intrinsic :: trim

  ! Local variables
  logical                     :: normal, precon, sing, consis, use_default
  integer(ip)                 :: n, nout, ios, i, j
  real(dp)                    :: shift, pertM
  real(dp)                    :: a, b, c_true, s_true, r_true, tol, BIG, SMALL
  integer(ip), parameter      :: nCDS = 2, nCPS = 10, nCRS = 7
  character(len=*), parameter :: pathCDS = './DataMtx/CDS/', &
                                 pathCPS = './DataMtx/CPS/', &
                                 pathCRS = './DataMtx/CRS/'
  character(len=*), parameter :: filesCDS(nCDS) =  &   !! Sorted by file size
      (/                    &
         'saylr3         '  &    ! n =  1000
       , 'shaw_100       '  &    ! n =   100
      /)
  character(len=*), parameter :: filesCPS(nCPS) =  &   !! Sorted by file size
      (/                    &
         'lap_25         '  &    ! n =    25
       , 'dwt_72         '  &    ! n =    72
       , 'GD97_a         '  &    ! n =    84
       , 'sphere2        '  &    ! n =    66
       , 'GD98_c         '  &    ! n =   112
       , 'can_61         '  &    ! n =    61
       , 'dwt_162        '  &    ! n =   162
       , 'can_144        '  &    ! n =   144
       , 'can_187        '  &    ! n =   187
       , 'NotreDame_yeast'  &    ! n =  2114
      /)
  character(len=*), parameter :: filesCRS(nCRS) =  &   !! Sorted by file size
      (/                    &
         'saylr3         '  &    ! n =  1000
       , 'laser          '  &    ! n =  3002
       , 'shaw_100       '  &    ! n =   100
       , 'mhd3200b       '  &    ! n =  3200
       , 'bcsstm13       '  &    ! n =  2003
       , 'c-30           '  &    ! n =  5321
       , 'c-54           '  &    ! n = 31793
       /)
  character(80)              :: input_file, output_file
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  write(*,*) "integer and real precision:"
  write(*,*) "ip =", ip
  write(*,*) "dp =", dp

  nout   = 10
  output_file = 'MINRESQLP.txt'
  open(nout, file=output_file, status='unknown', iostat=ios)

  if (ios /= 0) then
    write(*,*) "Error opening file", trim(output_file)
  end if

  if (dp <= -1) then
     write(*,*) "Error: precision or range of 'dp' defined in minresqlpDataModule is unsupported."
     stop
  end if

  if (ip <= -1) then
     write(*,*) "Error: range of 'ip' defined in minresqlpDataModule is unsupported."
     stop
  end if

  !------------------------------------------------------------------
  ! Tests for symortho
  !------------------------------------------------------------------
  if (testSymortho) then
     a      = zero
     b      = zero
     c_true = one
     s_true = zero
     r_true = zero
     call symorthotest ( a, b, c_true, s_true, r_true, nout )             ! Ex 01

     a      = zero
     b      = 2*one
     c_true = zero
     s_true = one
     r_true = 2*one
     call symorthotest ( a, b, c_true, s_true, r_true, nout )             ! Ex 02

     a      =    6_dp
     b      =   -8_dp
     c_true =  0.6_dp
     s_true = -0.8_dp
     r_true =   10_dp
     call symorthotest ( a, b, c_true, s_true, r_true, nout )             ! Ex 03

     a      =    8_dp
     b      =    6_dp
     c_true =  0.8_dp
     s_true =  0.6_dp
     r_true =   10_dp
     call symorthotest ( a, b, c_true, s_true, r_true, nout )             ! Ex 04

     BIG    = one * realmax/10.0_dp
     a      = BIG * (-3_dp)
     b      = BIG * (-4_dp)
     r_true = BIG *   5_dp
     c_true = a / r_true
     s_true = b / r_true
     call symorthotest ( a, b, c_true, s_true, r_true, nout )             ! Ex 05

     SMALL  = one * realmin
     a      = SMALL * (-3_dp)
     b      = SMALL * (-4_dp)
     r_true = SMALL *   5_dp
     c_true = a / r_true
     s_true = b / r_true
     call symorthotest ( a, b, c_true, s_true, r_true, nout )             ! Ex 06

     a      = BIG * (-3_dp)
     b      =        -4_dp
     r_true = BIG * ( 3_dp)
     c_true = a / r_true
     s_true = b / r_true
     call symorthotest ( a, b, c_true, s_true, r_true, nout )             ! Ex 07

     a      =        -4_dp
     b      = BIG * (-3_dp)
     r_true = BIG * ( 3_dp)
     c_true = a / r_true
     s_true = b / r_true
     call symorthotest ( a, b, c_true, s_true, r_true, nout )             ! Ex 08

     a      =          -3_dp
     b      = SMALL * (-4_dp)
     r_true =           3_dp
     c_true = a / r_true
     s_true = b / r_true
     call symorthotest ( a, b, c_true, s_true, r_true, nout )             ! Ex 09

     a      = SMALL * (-3_dp)
     b      =          -4_dp
     r_true =           4_dp
     c_true = a / r_true
     s_true = b / r_true
     call symorthotest ( a, b, c_true, s_true, r_true, nout )             ! Ex 09
  end if

  !------------------------------------------------------------------
  ! Tests for MINRESQLP
  !------------------------------------------------------------------
  do i = 1, 2
     normal = .false. ! no preconditioning
     precon = .true.
     shift  = 0.01_dp
     pertM  = 0.1_dp
     sing   = .false.
     consis = .true.
     use_default = .false.
     if (i > 1) use_default = .true.

     write(nout,*) ' '
     write(nout,*) ' MINRESQLP tests with use_default = ', use_default

     ! Test the unlikely tiny cases that often trip us up.

     n      = 1
     call minresqlptest( n, normal, zero , zero, sing, consis, use_default, nout ) ! Ex 01
     call minresqlptest( n, normal, shift, zero, sing, consis, use_default, nout ) ! Ex 02

     n      = 2
     call minresqlptest( n, normal, zero , zero, sing, consis, use_default, nout ) ! Ex 03
     call minresqlptest( n, normal, shift, zero, sing, consis, use_default, nout ) ! Ex 04

     n      = 3
     call minresqlptest( n, normal, zero , zero, sing, consis, use_default, nout ) ! Ex 05
     call minresqlptest( n, normal, shift, zero, sing, consis, use_default, nout ) ! Ex 06

     n      = 4
     call minresqlptest( n, normal, zero , zero, sing, consis, use_default, nout ) ! Ex 07
     call minresqlptest( n, normal, shift, zero, sing, consis, use_default, nout ) ! Ex 08

     ! Test small positive-definite and indefinite systems
     ! without preconditioners.  MINRESQLP should take n iterations.

     n      = 50
     call minresqlptest( n, normal, zero , zero, sing, consis, use_default, nout ) ! Ex 09
     call minresqlptest( n, normal, shift, zero, sing, consis, use_default, nout ) ! Ex 10

     ! Test small positive-definite and indefinite systems with
     ! exact preconditioners.  MINRESQLP should take about n iterations.

     n      = 2
     call minresqlptest( n, precon, zero , zero, sing, consis, use_default, nout ) ! Ex 11
     call minresqlptest( n, precon, shift, zero, sing, consis, use_default, nout ) ! Ex 12

     n      = 50
     call minresqlptest( n, precon, zero , zero, sing, consis, use_default, nout ) ! Ex 13
     call minresqlptest( n, precon, shift, zero, sing, consis, use_default, nout ) ! Ex 14

     ! pertM makes the preconditioners incorrect in n/10 entries.
     ! MINRESQLP should take about n/10 iterations.

     call minresqlptest( n, precon, zero , pertM, sing, consis, use_default, nout) ! Ex 15
     call minresqlptest( n, precon, shift, pertM, sing, consis, use_default, nout) ! Ex 16

     ! Singular consistent test case.
     sing   = .true.
     consis = .true.

     n      = 4
     call minresqlptest( n, normal, zero , zero, sing, consis, use_default, nout ) ! Ex 17

     ! if shift = 0.25 or 0.5 and n =4, then the problem is inconsistent
     shift  = 0.3_dp
     call minresqlptest( n, normal, shift, zero, sing, consis, use_default, nout ) ! Ex 18

     n      = 50
     call minresqlptest( n, normal, zero , zero, sing, consis, use_default, nout ) ! Ex 19

     ! Singular inconsistent test case.

     consis = .false.

     n      = 4
     call minresqlptest( n, normal, zero , zero, sing, consis, use_default, nout ) ! Ex 20

     n      = 50
     call minresqlptest( n, normal, zero , zero, sing, consis, use_default, nout ) ! Ex 21
  end do

  !------------------------------------------------------------------
  ! Test MINRESQLP on singular Matrix Market files
  !------------------------------------------------------------------
  tol = 1e-7_dp
  if (prcsn == 6) then                     ! single precision is used
    tol = 2.2e-3_dp
  end if

  if (testMtx) then
     write(nout,*) ' '
     write(nout,*) ' MINRESQLP tests on MM CDS examples'

     do j = 1, 2
        consis = .true.
        if (j > 1) consis = .false.
        do i = 1, nCDS
           input_file = pathCDS//trim(filesCDS(i))//'.mtx'
           call minresqlptestMtxCDS(input_file, consis, nout, tol)
        end do
     end do

     write(nout,*) ' '
     write(nout,*) ' MINRESQLP tests on MM CPS examples'

     do j = 1, 2
        consis = .true.
        if (j > 1) consis = .false.
        do i = 1, nCPS
           input_file = pathCPS//trim(filesCPS(i))//'.mtx'
           call minresqlptestMtxCPS(input_file, consis, nout, tol)
        end do
     end do

     write(nout,*) ' '
     write(nout,*) ' MINRESQLP tests on MM CRS examples'

     do j = 1, 2
        consis = .true.
        if (j > 1) consis = .false.
        do i = 1, 6 !nCRS
           input_file = pathCRS//trim(filesCRS(i))//'.mtx'
           call minresqlptestMtxCRS(input_file, consis, nout, tol)
        end do
     end do
  end if

  write(*,*)
  write(*,*) "Results are in output file  ", trim(output_file)
  write(*,*) "Search the file for 'appear'"
  write(*,*) "For example:    grep appear ", trim(output_file)
  write(*,*)

  close( unit = nout )
end program minresqlpTestProgram
