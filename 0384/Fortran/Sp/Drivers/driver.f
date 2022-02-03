      program main

c*******************************************************************************
c
cc TOMS384_PRB tests the TOMS384 routine SYMQR.
c
c  Modified:
c
c    09 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS384_PRB'
      write ( *, '(a)' ) '  Test SYMQR, which computes eigenvalues'
      write ( *, '(a)' ) '  and eigenvectors of a symmetric matrix.'

      call test01

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS384_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      stop
      end
      subroutine test01

c*******************************************************************************
c
cc TEST01 tests SYMQR on a full matrix.
c
c  Modified:
c
c    09 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 5 )

      real a(n,n)
      logical abscnv
      real angle
      real d(n)
      real e(n)
      real eps
      real exact
      integer fail
      integer i
      integer j
      real k0
      integer na
      real pi
      logical trd
      logical vec

      pi = 3.14159265E+00

      k0 = 0.0E+00
      na = n
      eps = 0.00001E+00
      abscnv = .false.
      vec = .true.
      trd = .false.

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
c
c  Set the values of the matrix A.
c
      do i = 1, n
        do j = 1, n
          a(i,j) = min ( i, j )
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The matrix A:'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,5f10.4)' ) ( a(i,j), j = 1, n )
      end do
c
c  Set the values of the right hand side vector B.
c
      call symqr ( a, d, e, k0, n, na, eps, abscnv, vec, trd,
     &  fail )

      if ( fail .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  FAIL was returned nonzero.'
        write ( *, '(a,i6)' ) '  FAIL = ', fail
        write ( *, '(a)' ) '  Only eigendata FAIL+1, ..., N can be'
        write ( *, '(a)' ) '  relied on.'
      end if
c
c  Print the results.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Eigenvalues:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    Computed          Exact'
      write ( *, '(a)' ) ' '

      do i = 1, n
        angle = real ( 2 * i - 1 ) * pi / real ( 2 * n + 1 )
        exact = 0.5E+00 / ( 1.0E+00 - cos ( angle ) )
        write ( *, '(2x,g14.6,2x,g14.6)' ) d(i), exact
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Eigenvector matrix:'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,5f10.4)' ) ( a(i,j), j = 1, n )
      end do

      return
      end
