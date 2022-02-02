! -- A Fully Portable High Performance Minimal --
! -- Storage Hybrid Format Cholesky Algorithm  --
!    F.G. Gustavson - IBM, USA;
!    J.K. Reid - RAL, UK;
!    J. Wasniewski - DTU, Denmark
!    March 20, 2006.

program lower_hybrid_speed

! For various values of n, nb, and mb, time the factorization to blocked hybrid 
! format and the solution of sets of equations using the factor.

  use block_hybrid_Cholesky
  implicit none
  character, parameter :: uplo = 'u' 
  integer, parameter ::  wp = kind(1.0d0)
  integer :: i, info, j, k, kmb, kn, knb, mb, nrhs, n, nb
  real(wp), allocatable :: a(:), ap(:), b(:,:), x(:,:)
  integer :: nc(11)=(/40,64,100,160,250,400,640,1000,1600,2500,4000/)
  integer :: nbc(4)=(/40,72,100,200/)
  integer :: mbc(3)=(/48,100,200/)
  real(wp) :: t, t1, t2, flops(0:3)  
  real(wp), parameter :: one = 1.0d0

  write(*,'(a)')        '              -----Mflops, uplo='//uplo//'------'
  write(*,'(a)')        '     n    nb  fact  -------solve-------'
  write(*,'(a,i3,2i6)') '                    mb =', mbc 

  do kn = 1,11
     n = nc(kn)
     nrhs = max(100,n/10)

! Allocate the arrays
     allocate(a(n*(n+1)/2),ap(n*(n+1)/2), b(n,nrhs), x(n,nrhs))

     do knb = 1,4
        nb = nbc(knb)
! Generate the lower-triangular matrix in the lower packed format
        call random_number(a)
        k = 1
        if (uplo=='u' .or. uplo=='U') then
           do i = 1,n
              a(k) = n
              k = k + i + 1
           end do
        else
           do i = n,1,-1
              a(k) = n
              k = k + i  
           end do
        end if
        call random_number(b)

! Transform the matrix to lower blocked hybrid format
        call PpHpp(uplo, n, a, info, nb )
        if (info /= 0) call terminate("PpHpp")

! Factorize the matrix
        t = 0
        do i = 1, huge(1)
          ap = a
          call cpu_time(t1)
          call HppTf( uplo, n, ap, info, nb )
          if (info /= 0) call terminate("HppTf")
          call cpu_time(t2)
          t = t + t2 - t1
          if (t>one) exit
        end do
        flops(0) = real(n)**3*i/(t*3d6)

! Solve the equations
        do kmb = 1,3
           mb = mbc(kmb)
           t = 0
           do i = 1, huge(1)
              x = b
              call cpu_time(t1)
              call HppTs( uplo, n, nrhs, ap, x, n, info, nb, mb )
              if (info /= 0) call terminate("HppTs1")
              call cpu_time(t2)
              t = t + t2 - t1
              if (t>one) exit
           end do
           flops(kmb) = real(n)**2*real(nrhs)*i*2d-6/t

! Check the residuals
           call HppPp(uplo, n, a, info, nb )
           do j = 1, nrhs
              call dspmv(uplo, n, one, a, x(:,j), 1, -one, b(:,j), 1)
              if ( maxval(abs(b(:,j))) > epsilon(one)*n ) then
                 write(*,*) ' Max. residual is ',  maxval(abs(b(:,j)))
                 stop
              end if
           end do

        end do
        write(*,'(2i6,f6.0,f10.0,2f6.0)') n, nb, flops

     end do

! Dellocate the arrays
     deallocate(a, ap, b, x)

  end do
      
contains

  subroutine terminate(name)
    character(*) name
    write(*,*) "Stopping after failure in ",name," with info = ", info
    stop
  end subroutine terminate

end program lower_hybrid_speed
