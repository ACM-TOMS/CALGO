      program main

c***********************************************************************
c
cc TOMS343_PRB tests TOMS343.
c
c  Modified:
c
c    20 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none


      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS343_PRB'
      write ( *, '(a)' ) '  Test TOMS algorithm 343, to compute'
      write ( *, '(a)' ) '  the eigenvalues and eigenvectors of'
      write ( *, '(a)' ) '  a real general matrix.'

      call test01

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS343_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '

      stop
      end
      subroutine test01

c***********************************************************************
c
cc TEST01 tests EIGENP.
c
c  Modified:
c
c    20 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer nm
      integer n1
      integer n2
      integer n3

      parameter ( nm = 5 )
      parameter ( n1 = 5 )
      parameter ( n2 = 4 )
      parameter ( n3 = 3 )

      real a(nm,nm)
      real a1(5,5)
      real a2(4,4)
      real a3(3,3)
      real evi(nm)
      real evr(nm)
      integer i
      integer indic(nm)
      integer j
      integer k
      integer n
      real t
      real veci(nm,nm)
      real vecr(nm,nm)

      data a1 /
     &  -0.5, 1.0, 0.0, 0.0, 0.0,
     &  -1.0, 0.0, 1.0, 0.0, 0.0,
     &  -1.0, 0.0, 0.0, 1.0, 0.0,
     &  -0.5, 0.0, 0.0, 0.0, 1.0,
     &  -1.0, 0.0, 0.0, 0.0, 0.0 /

      data a2 /
     &  -2.0, -7.0,  0.0, -1.0,
     &   1.0, -5.0, -1.0,  0.0,
     &   1.0, -2.0, -3.0, -1.0,
     &   1.0, -4.0, -2.0,  0.0 /

      data a3 /
     &   1.00, 0.10, 0.00,
     &   0.00, 1.00, 1.00,
     &   0.01, 0.00, 1.00 /

      t = 24.0E+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  Test EIGENP, which computes eigenvalues and'
      write ( *, '(a)' ) '  eigenvectors of a general real matrix.'
      write ( *, '(a)' ) ' '

      do k = 1, 3

        if ( k .eq. 1 ) then

          n = n1

          do i = 1, n
            do j = 1, n
              a(i,j) = a1(i,j)
            end do
          end do
      
        else if ( k .eq. 2 ) then

          n = n2

          do i = 1, n
            do j = 1, n
              a(i,j) = a2(i,j)
            end do
          end do

        else if ( k .eq. 3 ) then

          n = n3

          do i = 1, n
            do j = 1, n
              a(i,j) = a3(i,j)
            end do
          end do
      
        end if

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Matrix A:'
        write ( *, '(a)' ) ' '

        do i = 1, n
          write ( *, '(5f10.4)' ) ( a(i,j), j = 1, n )
        end do

        call eigenp ( n, nm, a, t, evr, evi, vecr, veci, indic )

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '       I   INDIC   Real part  Imag part'
        write ( *, '(a)' ) ' '

        do i = 1, n
          write ( *, '(2x,i6,2x,i6,2x,g14.6,2x,g14.6)' ) 
     &      i, indic(i), evr(i), evi(i)
        end do

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Real parts of eigenvectors:'
        write ( *, '(a)' ) ' '

        do i = 1, n
          write ( *, '(5f10.4)' ) ( vecr(i,j), j = 1, n )
        end do

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Imaginary parts of eigenvectors:'
        write ( *, '(a)' ) ' '

        do i = 1, n
          write ( *, '(5f10.4)' ) ( veci(i,j), j = 1, n )
        end do

      end do

      return
      end
