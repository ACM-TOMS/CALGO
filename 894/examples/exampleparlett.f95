! An example for the modified block Schur--Parlett method

program exampleparlett
  use scalesquare,   only : sasqtrphif
  use blockdecomp,   only : bspqtr1dcmpf, bspblock
  use schurparlett,  only : bspqtrphif

  implicit none

  integer,          parameter :: n = 128, upto = 5
  double precision, parameter :: one = 1
  double precision :: a(n,n), u(n,n,0:upto), v(n,n,0:upto)
  double precision :: err(2), nrm
  type(bspblock)   :: info
  integer          :: j

  print '(a)', '# This program compares the block Schur--Parlett algorithm with'
  print '(a)', '# the scaling and squaring method. Relative difference will be'
  print '(a)', '# printed.'
  print '(a)', '# phi_j   difference'

  a = generate_matrix(n)

  info = bspqtr1dcmpf(a)
  u = bspqtrphif(upto, one, a, info)    ! Schur--Parlett algorithm
  v = sasqtrphif(upto, one, a)          ! Scaling and squaring method

  ! comparison

  do j=0,upto
    nrm = norm1(v(:,:,j))
    err(1) = norm1(u(:,:,j)-v(:,:,j))/nrm
    print '(i5,e17.7)', j, err(1)
  end do

contains

  ! Generateing an quasi-upper triangular matrix.
  function generate_matrix(n)
    integer, intent(in) :: n
    double precision    :: radius, center
    double precision    :: generate_matrix(n,n), temp
    integer             :: i, j, subd(n)

    ! off-diagonal entries

    generate_matrix = 0
    do i=1,n
      do j=i+1,n
        call random_number(temp)
        generate_matrix(i,j) = (temp*2-1)*2
      end do
    end do

    ! a sub-diagonal pattern

    subd = 0
    call random_number(temp)
    if (one/2 < temp) subd(1) = 1
    do i=2,n-1
      if (subd(i-1) == 0) then
        call random_number(temp)
        if (one/2 < temp) subd(i) = 1
      end if
    end do

    ! eigenvalues

    radius = one * 20; center = one * 2
    i = 1
    do while (i <= n)
      generate_matrix(i,i) = - (i-1) * one / (one * n) * radius + center
      if (subd(i) == 1) then
        call random_number(temp)                        ! a conjugate pair
        generate_matrix(i  ,i+1) =  (temp*2-1)*2
        generate_matrix(i+1,i  ) = -generate_matrix(i  ,i+1)
        generate_matrix(i+1,i+1) =  generate_matrix(i  ,i  )
        i = i + 2
      else
        i = i + 1
      end if
    end do
  end function generate_matrix

  ! the matrix 1-norm
  pure function norm1(a)
    double precision, intent(in) :: a(:,:)
    double precision :: norm1, v(size(a,1))
    integer          :: i
    do i=1,size(a,1); v(i) = sum(abs(a(:,i))); end do
    norm1 = maxval(v)
  end function norm1

end program exampleparlett

