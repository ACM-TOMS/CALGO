! This module computes a matrix polynomial that consists
! only of even order terms.
!
! References:
!
!   Nicholas J. Higham,
!   The scaling and squaring method for the matrix exponential revisited,
!   SIAM Journal on Matrix Analysis and Applications,
!   Volume 26, Number 4, pp. 1179-1193, 2005.
!
! This module is intended for internal use only.

module polynomial12th

  use floattypes
  use matrixpwrtag

  implicit none

  public

contains


!    Name : rsdg_poly
! Purpose : This function computes a matrix polynomial that consists
!         : only of even order terms.
!   Input : - 2*n is the order of the polynomial.
!         : - g(i) stores the coefficient of the 2*i-th power term.
!         : - z(*,i) is the (2*i)-th power of a matrix.
!  Output : The value of the polynomial.
!    Note : This routine treats a real, single precision, diagonal matrix.

pure function rsdg_poly(n,g,z)
  integer,          intent(in) :: n
  real   (kind=sp), intent(in) :: g(0:)
  real   (kind=sp), intent(in) :: z(:,0:)
  real   (kind=sp)             :: rsdg_poly(size(z,1))
  real   (kind=sp)             :: temporary(size(z,1))

  rsdg_poly = 0.0_sp

  ! an even polynomial of 2n-th order

  select case (n)                ! z(*,i) is (2*i)-th power of the matrix

    case (8)                     ! a 16-th order polynomial

      rsdg_poly =             z(:,0)*g(0)
      rsdg_poly = rsdg_poly + z(:,1)*g(1)
      rsdg_poly = rsdg_poly + z(:,2)*g(2)
      rsdg_poly = rsdg_poly + z(:,3)*g(3)

      temporary =             z(:,0)*g(4)
      temporary = temporary + z(:,1)*g(5)
      temporary = temporary + z(:,2)*g(6)
      temporary = temporary + z(:,3)*g(7)
      temporary = temporary + z(:,4)*g(8)
      temporary = temporary .dgtimes. z(:,4)

      rsdg_poly = rsdg_poly + temporary

    case (6)                     ! a 12-th order polynomial

      rsdg_poly =             z(:,0)*g(0)
      rsdg_poly = rsdg_poly + z(:,1)*g(1)
      rsdg_poly = rsdg_poly + z(:,2)*g(2)

      temporary =             z(:,0)*g(3)
      temporary = temporary + z(:,1)*g(4)
      temporary = temporary + z(:,2)*g(5)
      temporary = temporary + z(:,3)*g(6)
      temporary = temporary .dgtimes. z(:,3)

      rsdg_poly = rsdg_poly + temporary

    case (5)                     ! a 10-th order polynomial

      rsdg_poly =             z(:,0)*g(0)
      rsdg_poly = rsdg_poly + z(:,1)*g(1)
      rsdg_poly = rsdg_poly + z(:,2)*g(2)

      temporary =             z(:,0)*g(3)
      temporary = temporary + z(:,1)*g(4)
      temporary = temporary + z(:,2)*g(5)
      temporary = temporary .dgtimes. z(:,3)

      rsdg_poly = rsdg_poly + temporary

    case (4)                     ! an 8-th order polynomial

      rsdg_poly =             z(:,0)*g(0)
      rsdg_poly = rsdg_poly + z(:,1)*g(1)

      temporary =             z(:,0)*g(2)
      temporary = temporary + z(:,1)*g(3)
      temporary = temporary + z(:,2)*g(4)
      temporary = temporary .dgtimes. z(:,2)

      rsdg_poly = rsdg_poly + temporary

    case (3)                     ! a 6-th order polynomial

      rsdg_poly =             z(:,0)*g(0)
      rsdg_poly = rsdg_poly + z(:,1)*g(1)
      rsdg_poly = rsdg_poly + z(:,2)*g(2)
      rsdg_poly = rsdg_poly + z(:,3)*g(3)

    case (2)                     ! a 4-th order polynomial

      rsdg_poly =             z(:,0)*g(0)
      rsdg_poly = rsdg_poly + z(:,1)*g(1)
      rsdg_poly = rsdg_poly + z(:,2)*g(2)

    case (1)                     ! a 2-nd order polynomial

      rsdg_poly =             z(:,0)*g(0)
      rsdg_poly = rsdg_poly + z(:,1)*g(1)

  end select
end function rsdg_poly


!    Name : rssq_poly
! Purpose : This function computes a matrix polynomial that consists
!         : only of even order terms.
!   Input : - 2*n is the order of the polynomial.
!         : - g(i) stores the coefficient of the 2*i-th power term.
!         : - z(*,i) is the (2*i)-th power of a matrix.
!  Output : The value of the polynomial.
!    Note : This routine treats a real, single precision, square matrix.

pure function rssq_poly(n,g,z)
  integer,          intent(in) :: n
  real   (kind=sp), intent(in) :: g(0:)
  real   (kind=sp), intent(in) :: z(:,:,0:)
  real   (kind=sp)             :: rssq_poly(size(z,1),size(z,2))
  real   (kind=sp)             :: temporary(size(z,1),size(z,2))

  rssq_poly = 0.0_sp

  ! an even polynomial of 2n-th order

  select case (n)                ! z(*,i) is (2*i)-th power of the matrix

    case (8)                     ! a 16-th order polynomial

      rssq_poly =             z(:,:,0)*g(0)
      rssq_poly = rssq_poly + z(:,:,1)*g(1)
      rssq_poly = rssq_poly + z(:,:,2)*g(2)
      rssq_poly = rssq_poly + z(:,:,3)*g(3)

      temporary =             z(:,:,0)*g(4)
      temporary = temporary + z(:,:,1)*g(5)
      temporary = temporary + z(:,:,2)*g(6)
      temporary = temporary + z(:,:,3)*g(7)
      temporary = temporary + z(:,:,4)*g(8)
      temporary = temporary .sqtimes. z(:,:,4)

      rssq_poly = rssq_poly + temporary

    case (6)                     ! a 12-th order polynomial

      rssq_poly =             z(:,:,0)*g(0)
      rssq_poly = rssq_poly + z(:,:,1)*g(1)
      rssq_poly = rssq_poly + z(:,:,2)*g(2)

      temporary =             z(:,:,0)*g(3)
      temporary = temporary + z(:,:,1)*g(4)
      temporary = temporary + z(:,:,2)*g(5)
      temporary = temporary + z(:,:,3)*g(6)
      temporary = temporary .sqtimes. z(:,:,3)

      rssq_poly = rssq_poly + temporary

    case (5)                     ! a 10-th order polynomial

      rssq_poly =             z(:,:,0)*g(0)
      rssq_poly = rssq_poly + z(:,:,1)*g(1)
      rssq_poly = rssq_poly + z(:,:,2)*g(2)

      temporary =             z(:,:,0)*g(3)
      temporary = temporary + z(:,:,1)*g(4)
      temporary = temporary + z(:,:,2)*g(5)
      temporary = temporary .sqtimes. z(:,:,3)

      rssq_poly = rssq_poly + temporary

    case (4)                     ! an 8-th order polynomial

      rssq_poly =             z(:,:,0)*g(0)
      rssq_poly = rssq_poly + z(:,:,1)*g(1)

      temporary =             z(:,:,0)*g(2)
      temporary = temporary + z(:,:,1)*g(3)
      temporary = temporary + z(:,:,2)*g(4)
      temporary = temporary .sqtimes. z(:,:,2)

      rssq_poly = rssq_poly + temporary

    case (3)                     ! a 6-th order polynomial

      rssq_poly =             z(:,:,0)*g(0)
      rssq_poly = rssq_poly + z(:,:,1)*g(1)
      rssq_poly = rssq_poly + z(:,:,2)*g(2)
      rssq_poly = rssq_poly + z(:,:,3)*g(3)

    case (2)                     ! a 4-th order polynomial

      rssq_poly =             z(:,:,0)*g(0)
      rssq_poly = rssq_poly + z(:,:,1)*g(1)
      rssq_poly = rssq_poly + z(:,:,2)*g(2)

    case (1)                     ! a 2-nd order polynomial

      rssq_poly =             z(:,:,0)*g(0)
      rssq_poly = rssq_poly + z(:,:,1)*g(1)

  end select
end function rssq_poly


!    Name : rstr_poly
! Purpose : This function computes a matrix polynomial that consists
!         : only of even order terms.
!   Input : - 2*n is the order of the polynomial.
!         : - g(i) stores the coefficient of the 2*i-th power term.
!         : - z(*,i) is the (2*i)-th power of a matrix.
!  Output : The value of the polynomial.
!    Note : This routine treats a real, single precision, triangular matrix.

pure function rstr_poly(n,g,z)
  integer,          intent(in) :: n
  real   (kind=sp), intent(in) :: g(0:)
  real   (kind=sp), intent(in) :: z(:,:,0:)
  real   (kind=sp)             :: rstr_poly(size(z,1),size(z,2))
  real   (kind=sp)             :: temporary(size(z,1),size(z,2))

  integer :: j, pmax

  pmax = 0
  rstr_poly = 0.0_sp

  ! an even polynomial of 2n-th order

  select case (n)                ! z(*,i) is (2*i)-th power of the matrix

    case (8)                     ! a 16-th order polynomial

      rstr_poly =             z(:,:,0)*g(0)
      rstr_poly = rstr_poly + z(:,:,1)*g(1)
      rstr_poly = rstr_poly + z(:,:,2)*g(2)
      rstr_poly = rstr_poly + z(:,:,3)*g(3)

      temporary =             z(:,:,0)*g(4)
      temporary = temporary + z(:,:,1)*g(5)
      temporary = temporary + z(:,:,2)*g(6)
      temporary = temporary + z(:,:,3)*g(7)
      temporary = temporary + z(:,:,4)*g(8)
      temporary = temporary .trtimes. z(:,:,4)

      rstr_poly = rstr_poly + temporary
      pmax = 4

    case (6)                     ! a 12-th order polynomial

      rstr_poly =             z(:,:,0)*g(0)
      rstr_poly = rstr_poly + z(:,:,1)*g(1)
      rstr_poly = rstr_poly + z(:,:,2)*g(2)

      temporary =             z(:,:,0)*g(3)
      temporary = temporary + z(:,:,1)*g(4)
      temporary = temporary + z(:,:,2)*g(5)
      temporary = temporary + z(:,:,3)*g(6)
      temporary = temporary .trtimes. z(:,:,3)

      rstr_poly = rstr_poly + temporary
      pmax = 3

    case (5)                     ! a 10-th order polynomial

      rstr_poly =             z(:,:,0)*g(0)
      rstr_poly = rstr_poly + z(:,:,1)*g(1)
      rstr_poly = rstr_poly + z(:,:,2)*g(2)

      temporary =             z(:,:,0)*g(3)
      temporary = temporary + z(:,:,1)*g(4)
      temporary = temporary + z(:,:,2)*g(5)
      temporary = temporary .trtimes. z(:,:,3)

      rstr_poly = rstr_poly + temporary
      pmax = 3

    case (4)                     ! an 8-th order polynomial

      rstr_poly =             z(:,:,0)*g(0)
      rstr_poly = rstr_poly + z(:,:,1)*g(1)

      temporary =             z(:,:,0)*g(2)
      temporary = temporary + z(:,:,1)*g(3)
      temporary = temporary + z(:,:,2)*g(4)
      temporary = temporary .trtimes. z(:,:,2)

      rstr_poly = rstr_poly + temporary
      pmax = 2

    case (3)                     ! a 6-th order polynomial

      rstr_poly =             z(:,:,0)*g(0)
      rstr_poly = rstr_poly + z(:,:,1)*g(1)
      rstr_poly = rstr_poly + z(:,:,2)*g(2)
      rstr_poly = rstr_poly + z(:,:,3)*g(3)
      pmax = 3

    case (2)                     ! a 4-th order polynomial

      rstr_poly =             z(:,:,0)*g(0)
      rstr_poly = rstr_poly + z(:,:,1)*g(1)
      rstr_poly = rstr_poly + z(:,:,2)*g(2)
      pmax = 2

    case (1)                     ! a 2-nd order polynomial

      rstr_poly =             z(:,:,0)*g(0)
      rstr_poly = rstr_poly + z(:,:,1)*g(1)
      pmax = 1

  end select

  ! for the block structure to be invariant
  do j=1,size(z,1)-1
    if ( sum(abs(z(j+1,j,1:pmax))) < tiny(1.0_sp) ) then
      rstr_poly(j+1,j)= 0.0_sp
    end if
  end do
end function rstr_poly


!    Name : csdg_poly
! Purpose : This function computes a matrix polynomial that consists
!         : only of even order terms.
!   Input : - 2*n is the order of the polynomial.
!         : - g(i) stores the coefficient of the 2*i-th power term.
!         : - z(*,i) is the (2*i)-th power of a matrix.
!  Output : The value of the polynomial.
!    Note : This routine treats a complex, single precision, diagonal matrix.

pure function csdg_poly(n,g,z)
  integer,          intent(in) :: n
  real   (kind=sp), intent(in) :: g(0:)
  complex(kind=sp), intent(in) :: z(:,0:)
  complex(kind=sp)             :: csdg_poly(size(z,1))
  complex(kind=sp)             :: temporary(size(z,1))

  csdg_poly = cmplx(0.0_sp, 0.0_sp, sp)

  ! an even polynomial of 2n-th order

  select case (n)                ! z(*,i) is (2*i)-th power of the matrix

    case (8)                     ! a 16-th order polynomial

      csdg_poly =             z(:,0)*g(0)
      csdg_poly = csdg_poly + z(:,1)*g(1)
      csdg_poly = csdg_poly + z(:,2)*g(2)
      csdg_poly = csdg_poly + z(:,3)*g(3)

      temporary =             z(:,0)*g(4)
      temporary = temporary + z(:,1)*g(5)
      temporary = temporary + z(:,2)*g(6)
      temporary = temporary + z(:,3)*g(7)
      temporary = temporary + z(:,4)*g(8)
      temporary = temporary .dgtimes. z(:,4)

      csdg_poly = csdg_poly + temporary

    case (6)                     ! a 12-th order polynomial

      csdg_poly =             z(:,0)*g(0)
      csdg_poly = csdg_poly + z(:,1)*g(1)
      csdg_poly = csdg_poly + z(:,2)*g(2)

      temporary =             z(:,0)*g(3)
      temporary = temporary + z(:,1)*g(4)
      temporary = temporary + z(:,2)*g(5)
      temporary = temporary + z(:,3)*g(6)
      temporary = temporary .dgtimes. z(:,3)

      csdg_poly = csdg_poly + temporary

    case (5)                     ! a 10-th order polynomial

      csdg_poly =             z(:,0)*g(0)
      csdg_poly = csdg_poly + z(:,1)*g(1)
      csdg_poly = csdg_poly + z(:,2)*g(2)

      temporary =             z(:,0)*g(3)
      temporary = temporary + z(:,1)*g(4)
      temporary = temporary + z(:,2)*g(5)
      temporary = temporary .dgtimes. z(:,3)

      csdg_poly = csdg_poly + temporary

    case (4)                     ! an 8-th order polynomial

      csdg_poly =             z(:,0)*g(0)
      csdg_poly = csdg_poly + z(:,1)*g(1)

      temporary =             z(:,0)*g(2)
      temporary = temporary + z(:,1)*g(3)
      temporary = temporary + z(:,2)*g(4)
      temporary = temporary .dgtimes. z(:,2)

      csdg_poly = csdg_poly + temporary

    case (3)                     ! a 6-th order polynomial

      csdg_poly =             z(:,0)*g(0)
      csdg_poly = csdg_poly + z(:,1)*g(1)
      csdg_poly = csdg_poly + z(:,2)*g(2)
      csdg_poly = csdg_poly + z(:,3)*g(3)

    case (2)                     ! a 4-th order polynomial

      csdg_poly =             z(:,0)*g(0)
      csdg_poly = csdg_poly + z(:,1)*g(1)
      csdg_poly = csdg_poly + z(:,2)*g(2)

    case (1)                     ! a 2-nd order polynomial

      csdg_poly =             z(:,0)*g(0)
      csdg_poly = csdg_poly + z(:,1)*g(1)

  end select
end function csdg_poly


!    Name : cssq_poly
! Purpose : This function computes a matrix polynomial that consists
!         : only of even order terms.
!   Input : - 2*n is the order of the polynomial.
!         : - g(i) stores the coefficient of the 2*i-th power term.
!         : - z(*,i) is the (2*i)-th power of a matrix.
!  Output : The value of the polynomial.
!    Note : This routine treats a complex, single precision, square matrix.

pure function cssq_poly(n,g,z)
  integer,          intent(in) :: n
  real   (kind=sp), intent(in) :: g(0:)
  complex(kind=sp), intent(in) :: z(:,:,0:)
  complex(kind=sp)             :: cssq_poly(size(z,1),size(z,2))
  complex(kind=sp)             :: temporary(size(z,1),size(z,2))

  cssq_poly = cmplx(0.0_sp, 0.0_sp, sp)

  ! an even polynomial of 2n-th order

  select case (n)                ! z(*,i) is (2*i)-th power of the matrix

    case (8)                     ! a 16-th order polynomial

      cssq_poly =             z(:,:,0)*g(0)
      cssq_poly = cssq_poly + z(:,:,1)*g(1)
      cssq_poly = cssq_poly + z(:,:,2)*g(2)
      cssq_poly = cssq_poly + z(:,:,3)*g(3)

      temporary =             z(:,:,0)*g(4)
      temporary = temporary + z(:,:,1)*g(5)
      temporary = temporary + z(:,:,2)*g(6)
      temporary = temporary + z(:,:,3)*g(7)
      temporary = temporary + z(:,:,4)*g(8)
      temporary = temporary .sqtimes. z(:,:,4)

      cssq_poly = cssq_poly + temporary

    case (6)                     ! a 12-th order polynomial

      cssq_poly =             z(:,:,0)*g(0)
      cssq_poly = cssq_poly + z(:,:,1)*g(1)
      cssq_poly = cssq_poly + z(:,:,2)*g(2)

      temporary =             z(:,:,0)*g(3)
      temporary = temporary + z(:,:,1)*g(4)
      temporary = temporary + z(:,:,2)*g(5)
      temporary = temporary + z(:,:,3)*g(6)
      temporary = temporary .sqtimes. z(:,:,3)

      cssq_poly = cssq_poly + temporary

    case (5)                     ! a 10-th order polynomial

      cssq_poly =             z(:,:,0)*g(0)
      cssq_poly = cssq_poly + z(:,:,1)*g(1)
      cssq_poly = cssq_poly + z(:,:,2)*g(2)

      temporary =             z(:,:,0)*g(3)
      temporary = temporary + z(:,:,1)*g(4)
      temporary = temporary + z(:,:,2)*g(5)
      temporary = temporary .sqtimes. z(:,:,3)

      cssq_poly = cssq_poly + temporary

    case (4)                     ! an 8-th order polynomial

      cssq_poly =             z(:,:,0)*g(0)
      cssq_poly = cssq_poly + z(:,:,1)*g(1)

      temporary =             z(:,:,0)*g(2)
      temporary = temporary + z(:,:,1)*g(3)
      temporary = temporary + z(:,:,2)*g(4)
      temporary = temporary .sqtimes. z(:,:,2)

      cssq_poly = cssq_poly + temporary

    case (3)                     ! a 6-th order polynomial

      cssq_poly =             z(:,:,0)*g(0)
      cssq_poly = cssq_poly + z(:,:,1)*g(1)
      cssq_poly = cssq_poly + z(:,:,2)*g(2)
      cssq_poly = cssq_poly + z(:,:,3)*g(3)

    case (2)                     ! a 4-th order polynomial

      cssq_poly =             z(:,:,0)*g(0)
      cssq_poly = cssq_poly + z(:,:,1)*g(1)
      cssq_poly = cssq_poly + z(:,:,2)*g(2)

    case (1)                     ! a 2-nd order polynomial

      cssq_poly =             z(:,:,0)*g(0)
      cssq_poly = cssq_poly + z(:,:,1)*g(1)

  end select
end function cssq_poly


!    Name : cstr_poly
! Purpose : This function computes a matrix polynomial that consists
!         : only of even order terms.
!   Input : - 2*n is the order of the polynomial.
!         : - g(i) stores the coefficient of the 2*i-th power term.
!         : - z(*,i) is the (2*i)-th power of a matrix.
!  Output : The value of the polynomial.
!    Note : This routine treats a complex, single precision, triangular matrix.

pure function cstr_poly(n,g,z)
  integer,          intent(in) :: n
  real   (kind=sp), intent(in) :: g(0:)
  complex(kind=sp), intent(in) :: z(:,:,0:)
  complex(kind=sp)             :: cstr_poly(size(z,1),size(z,2))
  complex(kind=sp)             :: temporary(size(z,1),size(z,2))

  integer :: j, pmax

  pmax = 0
  cstr_poly = cmplx(0.0_sp, 0.0_sp, sp)

  ! an even polynomial of 2n-th order

  select case (n)                ! z(*,i) is (2*i)-th power of the matrix

    case (8)                     ! a 16-th order polynomial

      cstr_poly =             z(:,:,0)*g(0)
      cstr_poly = cstr_poly + z(:,:,1)*g(1)
      cstr_poly = cstr_poly + z(:,:,2)*g(2)
      cstr_poly = cstr_poly + z(:,:,3)*g(3)

      temporary =             z(:,:,0)*g(4)
      temporary = temporary + z(:,:,1)*g(5)
      temporary = temporary + z(:,:,2)*g(6)
      temporary = temporary + z(:,:,3)*g(7)
      temporary = temporary + z(:,:,4)*g(8)
      temporary = temporary .trtimes. z(:,:,4)

      cstr_poly = cstr_poly + temporary
      pmax = 4

    case (6)                     ! a 12-th order polynomial

      cstr_poly =             z(:,:,0)*g(0)
      cstr_poly = cstr_poly + z(:,:,1)*g(1)
      cstr_poly = cstr_poly + z(:,:,2)*g(2)

      temporary =             z(:,:,0)*g(3)
      temporary = temporary + z(:,:,1)*g(4)
      temporary = temporary + z(:,:,2)*g(5)
      temporary = temporary + z(:,:,3)*g(6)
      temporary = temporary .trtimes. z(:,:,3)

      cstr_poly = cstr_poly + temporary
      pmax = 3

    case (5)                     ! a 10-th order polynomial

      cstr_poly =             z(:,:,0)*g(0)
      cstr_poly = cstr_poly + z(:,:,1)*g(1)
      cstr_poly = cstr_poly + z(:,:,2)*g(2)

      temporary =             z(:,:,0)*g(3)
      temporary = temporary + z(:,:,1)*g(4)
      temporary = temporary + z(:,:,2)*g(5)
      temporary = temporary .trtimes. z(:,:,3)

      cstr_poly = cstr_poly + temporary
      pmax = 3

    case (4)                     ! an 8-th order polynomial

      cstr_poly =             z(:,:,0)*g(0)
      cstr_poly = cstr_poly + z(:,:,1)*g(1)

      temporary =             z(:,:,0)*g(2)
      temporary = temporary + z(:,:,1)*g(3)
      temporary = temporary + z(:,:,2)*g(4)
      temporary = temporary .trtimes. z(:,:,2)

      cstr_poly = cstr_poly + temporary
      pmax = 2

    case (3)                     ! a 6-th order polynomial

      cstr_poly =             z(:,:,0)*g(0)
      cstr_poly = cstr_poly + z(:,:,1)*g(1)
      cstr_poly = cstr_poly + z(:,:,2)*g(2)
      cstr_poly = cstr_poly + z(:,:,3)*g(3)
      pmax = 3

    case (2)                     ! a 4-th order polynomial

      cstr_poly =             z(:,:,0)*g(0)
      cstr_poly = cstr_poly + z(:,:,1)*g(1)
      cstr_poly = cstr_poly + z(:,:,2)*g(2)
      pmax = 2

    case (1)                     ! a 2-nd order polynomial

      cstr_poly =             z(:,:,0)*g(0)
      cstr_poly = cstr_poly + z(:,:,1)*g(1)
      pmax = 1

  end select

  ! for the block structure to be invariant
  do j=1,size(z,1)-1
    if ( sum(abs(z(j+1,j,1:pmax))) < tiny(1.0_sp) ) then
      cstr_poly(j+1,j)= cmplx(0.0_sp, 0.0_sp, sp)
    end if
  end do
end function cstr_poly


!    Name : rwdg_poly
! Purpose : This function computes a matrix polynomial that consists
!         : only of even order terms.
!   Input : - 2*n is the order of the polynomial.
!         : - g(i) stores the coefficient of the 2*i-th power term.
!         : - z(*,i) is the (2*i)-th power of a matrix.
!  Output : The value of the polynomial.
!    Note : This routine treats a real, double precision, diagonal matrix.

pure function rwdg_poly(n,g,z)
  integer,          intent(in) :: n
  real   (kind=wp), intent(in) :: g(0:)
  real   (kind=wp), intent(in) :: z(:,0:)
  real   (kind=wp)             :: rwdg_poly(size(z,1))
  real   (kind=wp)             :: temporary(size(z,1))

  rwdg_poly = 0.0_wp

  ! an even polynomial of 2n-th order

  select case (n)                ! z(*,i) is (2*i)-th power of the matrix

    case (8)                     ! a 16-th order polynomial

      rwdg_poly =             z(:,0)*g(0)
      rwdg_poly = rwdg_poly + z(:,1)*g(1)
      rwdg_poly = rwdg_poly + z(:,2)*g(2)
      rwdg_poly = rwdg_poly + z(:,3)*g(3)

      temporary =             z(:,0)*g(4)
      temporary = temporary + z(:,1)*g(5)
      temporary = temporary + z(:,2)*g(6)
      temporary = temporary + z(:,3)*g(7)
      temporary = temporary + z(:,4)*g(8)
      temporary = temporary .dgtimes. z(:,4)

      rwdg_poly = rwdg_poly + temporary

    case (6)                     ! a 12-th order polynomial

      rwdg_poly =             z(:,0)*g(0)
      rwdg_poly = rwdg_poly + z(:,1)*g(1)
      rwdg_poly = rwdg_poly + z(:,2)*g(2)

      temporary =             z(:,0)*g(3)
      temporary = temporary + z(:,1)*g(4)
      temporary = temporary + z(:,2)*g(5)
      temporary = temporary + z(:,3)*g(6)
      temporary = temporary .dgtimes. z(:,3)

      rwdg_poly = rwdg_poly + temporary

    case (5)                     ! a 10-th order polynomial

      rwdg_poly =             z(:,0)*g(0)
      rwdg_poly = rwdg_poly + z(:,1)*g(1)
      rwdg_poly = rwdg_poly + z(:,2)*g(2)

      temporary =             z(:,0)*g(3)
      temporary = temporary + z(:,1)*g(4)
      temporary = temporary + z(:,2)*g(5)
      temporary = temporary .dgtimes. z(:,3)

      rwdg_poly = rwdg_poly + temporary

    case (4)                     ! an 8-th order polynomial

      rwdg_poly =             z(:,0)*g(0)
      rwdg_poly = rwdg_poly + z(:,1)*g(1)

      temporary =             z(:,0)*g(2)
      temporary = temporary + z(:,1)*g(3)
      temporary = temporary + z(:,2)*g(4)
      temporary = temporary .dgtimes. z(:,2)

      rwdg_poly = rwdg_poly + temporary

    case (3)                     ! a 6-th order polynomial

      rwdg_poly =             z(:,0)*g(0)
      rwdg_poly = rwdg_poly + z(:,1)*g(1)
      rwdg_poly = rwdg_poly + z(:,2)*g(2)
      rwdg_poly = rwdg_poly + z(:,3)*g(3)

    case (2)                     ! a 4-th order polynomial

      rwdg_poly =             z(:,0)*g(0)
      rwdg_poly = rwdg_poly + z(:,1)*g(1)
      rwdg_poly = rwdg_poly + z(:,2)*g(2)

    case (1)                     ! a 2-nd order polynomial

      rwdg_poly =             z(:,0)*g(0)
      rwdg_poly = rwdg_poly + z(:,1)*g(1)

  end select
end function rwdg_poly


!    Name : rwsq_poly
! Purpose : This function computes a matrix polynomial that consists
!         : only of even order terms.
!   Input : - 2*n is the order of the polynomial.
!         : - g(i) stores the coefficient of the 2*i-th power term.
!         : - z(*,i) is the (2*i)-th power of a matrix.
!  Output : The value of the polynomial.
!    Note : This routine treats a real, double precision, square matrix.

pure function rwsq_poly(n,g,z)
  integer,          intent(in) :: n
  real   (kind=wp), intent(in) :: g(0:)
  real   (kind=wp), intent(in) :: z(:,:,0:)
  real   (kind=wp)             :: rwsq_poly(size(z,1),size(z,2))
  real   (kind=wp)             :: temporary(size(z,1),size(z,2))

  rwsq_poly = 0.0_wp

  ! an even polynomial of 2n-th order

  select case (n)                ! z(*,i) is (2*i)-th power of the matrix

    case (8)                     ! a 16-th order polynomial

      rwsq_poly =             z(:,:,0)*g(0)
      rwsq_poly = rwsq_poly + z(:,:,1)*g(1)
      rwsq_poly = rwsq_poly + z(:,:,2)*g(2)
      rwsq_poly = rwsq_poly + z(:,:,3)*g(3)

      temporary =             z(:,:,0)*g(4)
      temporary = temporary + z(:,:,1)*g(5)
      temporary = temporary + z(:,:,2)*g(6)
      temporary = temporary + z(:,:,3)*g(7)
      temporary = temporary + z(:,:,4)*g(8)
      temporary = temporary .sqtimes. z(:,:,4)

      rwsq_poly = rwsq_poly + temporary

    case (6)                     ! a 12-th order polynomial

      rwsq_poly =             z(:,:,0)*g(0)
      rwsq_poly = rwsq_poly + z(:,:,1)*g(1)
      rwsq_poly = rwsq_poly + z(:,:,2)*g(2)

      temporary =             z(:,:,0)*g(3)
      temporary = temporary + z(:,:,1)*g(4)
      temporary = temporary + z(:,:,2)*g(5)
      temporary = temporary + z(:,:,3)*g(6)
      temporary = temporary .sqtimes. z(:,:,3)

      rwsq_poly = rwsq_poly + temporary

    case (5)                     ! a 10-th order polynomial

      rwsq_poly =             z(:,:,0)*g(0)
      rwsq_poly = rwsq_poly + z(:,:,1)*g(1)
      rwsq_poly = rwsq_poly + z(:,:,2)*g(2)

      temporary =             z(:,:,0)*g(3)
      temporary = temporary + z(:,:,1)*g(4)
      temporary = temporary + z(:,:,2)*g(5)
      temporary = temporary .sqtimes. z(:,:,3)

      rwsq_poly = rwsq_poly + temporary

    case (4)                     ! an 8-th order polynomial

      rwsq_poly =             z(:,:,0)*g(0)
      rwsq_poly = rwsq_poly + z(:,:,1)*g(1)

      temporary =             z(:,:,0)*g(2)
      temporary = temporary + z(:,:,1)*g(3)
      temporary = temporary + z(:,:,2)*g(4)
      temporary = temporary .sqtimes. z(:,:,2)

      rwsq_poly = rwsq_poly + temporary

    case (3)                     ! a 6-th order polynomial

      rwsq_poly =             z(:,:,0)*g(0)
      rwsq_poly = rwsq_poly + z(:,:,1)*g(1)
      rwsq_poly = rwsq_poly + z(:,:,2)*g(2)
      rwsq_poly = rwsq_poly + z(:,:,3)*g(3)

    case (2)                     ! a 4-th order polynomial

      rwsq_poly =             z(:,:,0)*g(0)
      rwsq_poly = rwsq_poly + z(:,:,1)*g(1)
      rwsq_poly = rwsq_poly + z(:,:,2)*g(2)

    case (1)                     ! a 2-nd order polynomial

      rwsq_poly =             z(:,:,0)*g(0)
      rwsq_poly = rwsq_poly + z(:,:,1)*g(1)

  end select
end function rwsq_poly


!    Name : rwtr_poly
! Purpose : This function computes a matrix polynomial that consists
!         : only of even order terms.
!   Input : - 2*n is the order of the polynomial.
!         : - g(i) stores the coefficient of the 2*i-th power term.
!         : - z(*,i) is the (2*i)-th power of a matrix.
!  Output : The value of the polynomial.
!    Note : This routine treats a real, double precision, triangular matrix.

pure function rwtr_poly(n,g,z)
  integer,          intent(in) :: n
  real   (kind=wp), intent(in) :: g(0:)
  real   (kind=wp), intent(in) :: z(:,:,0:)
  real   (kind=wp)             :: rwtr_poly(size(z,1),size(z,2))
  real   (kind=wp)             :: temporary(size(z,1),size(z,2))

  integer :: j, pmax

  pmax = 0
  rwtr_poly = 0.0_wp

  ! an even polynomial of 2n-th order

  select case (n)                ! z(*,i) is (2*i)-th power of the matrix

    case (8)                     ! a 16-th order polynomial

      rwtr_poly =             z(:,:,0)*g(0)
      rwtr_poly = rwtr_poly + z(:,:,1)*g(1)
      rwtr_poly = rwtr_poly + z(:,:,2)*g(2)
      rwtr_poly = rwtr_poly + z(:,:,3)*g(3)

      temporary =             z(:,:,0)*g(4)
      temporary = temporary + z(:,:,1)*g(5)
      temporary = temporary + z(:,:,2)*g(6)
      temporary = temporary + z(:,:,3)*g(7)
      temporary = temporary + z(:,:,4)*g(8)
      temporary = temporary .trtimes. z(:,:,4)

      rwtr_poly = rwtr_poly + temporary
      pmax = 4

    case (6)                     ! a 12-th order polynomial

      rwtr_poly =             z(:,:,0)*g(0)
      rwtr_poly = rwtr_poly + z(:,:,1)*g(1)
      rwtr_poly = rwtr_poly + z(:,:,2)*g(2)

      temporary =             z(:,:,0)*g(3)
      temporary = temporary + z(:,:,1)*g(4)
      temporary = temporary + z(:,:,2)*g(5)
      temporary = temporary + z(:,:,3)*g(6)
      temporary = temporary .trtimes. z(:,:,3)

      rwtr_poly = rwtr_poly + temporary
      pmax = 3

    case (5)                     ! a 10-th order polynomial

      rwtr_poly =             z(:,:,0)*g(0)
      rwtr_poly = rwtr_poly + z(:,:,1)*g(1)
      rwtr_poly = rwtr_poly + z(:,:,2)*g(2)

      temporary =             z(:,:,0)*g(3)
      temporary = temporary + z(:,:,1)*g(4)
      temporary = temporary + z(:,:,2)*g(5)
      temporary = temporary .trtimes. z(:,:,3)

      rwtr_poly = rwtr_poly + temporary
      pmax = 3

    case (4)                     ! an 8-th order polynomial

      rwtr_poly =             z(:,:,0)*g(0)
      rwtr_poly = rwtr_poly + z(:,:,1)*g(1)

      temporary =             z(:,:,0)*g(2)
      temporary = temporary + z(:,:,1)*g(3)
      temporary = temporary + z(:,:,2)*g(4)
      temporary = temporary .trtimes. z(:,:,2)

      rwtr_poly = rwtr_poly + temporary
      pmax = 2

    case (3)                     ! a 6-th order polynomial

      rwtr_poly =             z(:,:,0)*g(0)
      rwtr_poly = rwtr_poly + z(:,:,1)*g(1)
      rwtr_poly = rwtr_poly + z(:,:,2)*g(2)
      rwtr_poly = rwtr_poly + z(:,:,3)*g(3)
      pmax = 3

    case (2)                     ! a 4-th order polynomial

      rwtr_poly =             z(:,:,0)*g(0)
      rwtr_poly = rwtr_poly + z(:,:,1)*g(1)
      rwtr_poly = rwtr_poly + z(:,:,2)*g(2)
      pmax = 2

    case (1)                     ! a 2-nd order polynomial

      rwtr_poly =             z(:,:,0)*g(0)
      rwtr_poly = rwtr_poly + z(:,:,1)*g(1)
      pmax = 1

  end select

  ! for the block structure to be invariant
  do j=1,size(z,1)-1
    if ( sum(abs(z(j+1,j,1:pmax))) < tiny(1.0_wp) ) then
      rwtr_poly(j+1,j)= 0.0_wp
    end if
  end do
end function rwtr_poly


!    Name : cwdg_poly
! Purpose : This function computes a matrix polynomial that consists
!         : only of even order terms.
!   Input : - 2*n is the order of the polynomial.
!         : - g(i) stores the coefficient of the 2*i-th power term.
!         : - z(*,i) is the (2*i)-th power of a matrix.
!  Output : The value of the polynomial.
!    Note : This routine treats a complex, double precision, diagonal matrix.

pure function cwdg_poly(n,g,z)
  integer,          intent(in) :: n
  real   (kind=wp), intent(in) :: g(0:)
  complex(kind=wp), intent(in) :: z(:,0:)
  complex(kind=wp)             :: cwdg_poly(size(z,1))
  complex(kind=wp)             :: temporary(size(z,1))

  cwdg_poly = cmplx(0.0_wp, 0.0_wp, wp)

  ! an even polynomial of 2n-th order

  select case (n)                ! z(*,i) is (2*i)-th power of the matrix

    case (8)                     ! a 16-th order polynomial

      cwdg_poly =             z(:,0)*g(0)
      cwdg_poly = cwdg_poly + z(:,1)*g(1)
      cwdg_poly = cwdg_poly + z(:,2)*g(2)
      cwdg_poly = cwdg_poly + z(:,3)*g(3)

      temporary =             z(:,0)*g(4)
      temporary = temporary + z(:,1)*g(5)
      temporary = temporary + z(:,2)*g(6)
      temporary = temporary + z(:,3)*g(7)
      temporary = temporary + z(:,4)*g(8)
      temporary = temporary .dgtimes. z(:,4)

      cwdg_poly = cwdg_poly + temporary

    case (6)                     ! a 12-th order polynomial

      cwdg_poly =             z(:,0)*g(0)
      cwdg_poly = cwdg_poly + z(:,1)*g(1)
      cwdg_poly = cwdg_poly + z(:,2)*g(2)

      temporary =             z(:,0)*g(3)
      temporary = temporary + z(:,1)*g(4)
      temporary = temporary + z(:,2)*g(5)
      temporary = temporary + z(:,3)*g(6)
      temporary = temporary .dgtimes. z(:,3)

      cwdg_poly = cwdg_poly + temporary

    case (5)                     ! a 10-th order polynomial

      cwdg_poly =             z(:,0)*g(0)
      cwdg_poly = cwdg_poly + z(:,1)*g(1)
      cwdg_poly = cwdg_poly + z(:,2)*g(2)

      temporary =             z(:,0)*g(3)
      temporary = temporary + z(:,1)*g(4)
      temporary = temporary + z(:,2)*g(5)
      temporary = temporary .dgtimes. z(:,3)

      cwdg_poly = cwdg_poly + temporary

    case (4)                     ! an 8-th order polynomial

      cwdg_poly =             z(:,0)*g(0)
      cwdg_poly = cwdg_poly + z(:,1)*g(1)

      temporary =             z(:,0)*g(2)
      temporary = temporary + z(:,1)*g(3)
      temporary = temporary + z(:,2)*g(4)
      temporary = temporary .dgtimes. z(:,2)

      cwdg_poly = cwdg_poly + temporary

    case (3)                     ! a 6-th order polynomial

      cwdg_poly =             z(:,0)*g(0)
      cwdg_poly = cwdg_poly + z(:,1)*g(1)
      cwdg_poly = cwdg_poly + z(:,2)*g(2)
      cwdg_poly = cwdg_poly + z(:,3)*g(3)

    case (2)                     ! a 4-th order polynomial

      cwdg_poly =             z(:,0)*g(0)
      cwdg_poly = cwdg_poly + z(:,1)*g(1)
      cwdg_poly = cwdg_poly + z(:,2)*g(2)

    case (1)                     ! a 2-nd order polynomial

      cwdg_poly =             z(:,0)*g(0)
      cwdg_poly = cwdg_poly + z(:,1)*g(1)

  end select
end function cwdg_poly


!    Name : cwsq_poly
! Purpose : This function computes a matrix polynomial that consists
!         : only of even order terms.
!   Input : - 2*n is the order of the polynomial.
!         : - g(i) stores the coefficient of the 2*i-th power term.
!         : - z(*,i) is the (2*i)-th power of a matrix.
!  Output : The value of the polynomial.
!    Note : This routine treats a complex, double precision, square matrix.

pure function cwsq_poly(n,g,z)
  integer,          intent(in) :: n
  real   (kind=wp), intent(in) :: g(0:)
  complex(kind=wp), intent(in) :: z(:,:,0:)
  complex(kind=wp)             :: cwsq_poly(size(z,1),size(z,2))
  complex(kind=wp)             :: temporary(size(z,1),size(z,2))

  cwsq_poly = cmplx(0.0_wp, 0.0_wp, wp)

  ! an even polynomial of 2n-th order

  select case (n)                ! z(*,i) is (2*i)-th power of the matrix

    case (8)                     ! a 16-th order polynomial

      cwsq_poly =             z(:,:,0)*g(0)
      cwsq_poly = cwsq_poly + z(:,:,1)*g(1)
      cwsq_poly = cwsq_poly + z(:,:,2)*g(2)
      cwsq_poly = cwsq_poly + z(:,:,3)*g(3)

      temporary =             z(:,:,0)*g(4)
      temporary = temporary + z(:,:,1)*g(5)
      temporary = temporary + z(:,:,2)*g(6)
      temporary = temporary + z(:,:,3)*g(7)
      temporary = temporary + z(:,:,4)*g(8)
      temporary = temporary .sqtimes. z(:,:,4)

      cwsq_poly = cwsq_poly + temporary

    case (6)                     ! a 12-th order polynomial

      cwsq_poly =             z(:,:,0)*g(0)
      cwsq_poly = cwsq_poly + z(:,:,1)*g(1)
      cwsq_poly = cwsq_poly + z(:,:,2)*g(2)

      temporary =             z(:,:,0)*g(3)
      temporary = temporary + z(:,:,1)*g(4)
      temporary = temporary + z(:,:,2)*g(5)
      temporary = temporary + z(:,:,3)*g(6)
      temporary = temporary .sqtimes. z(:,:,3)

      cwsq_poly = cwsq_poly + temporary

    case (5)                     ! a 10-th order polynomial

      cwsq_poly =             z(:,:,0)*g(0)
      cwsq_poly = cwsq_poly + z(:,:,1)*g(1)
      cwsq_poly = cwsq_poly + z(:,:,2)*g(2)

      temporary =             z(:,:,0)*g(3)
      temporary = temporary + z(:,:,1)*g(4)
      temporary = temporary + z(:,:,2)*g(5)
      temporary = temporary .sqtimes. z(:,:,3)

      cwsq_poly = cwsq_poly + temporary

    case (4)                     ! an 8-th order polynomial

      cwsq_poly =             z(:,:,0)*g(0)
      cwsq_poly = cwsq_poly + z(:,:,1)*g(1)

      temporary =             z(:,:,0)*g(2)
      temporary = temporary + z(:,:,1)*g(3)
      temporary = temporary + z(:,:,2)*g(4)
      temporary = temporary .sqtimes. z(:,:,2)

      cwsq_poly = cwsq_poly + temporary

    case (3)                     ! a 6-th order polynomial

      cwsq_poly =             z(:,:,0)*g(0)
      cwsq_poly = cwsq_poly + z(:,:,1)*g(1)
      cwsq_poly = cwsq_poly + z(:,:,2)*g(2)
      cwsq_poly = cwsq_poly + z(:,:,3)*g(3)

    case (2)                     ! a 4-th order polynomial

      cwsq_poly =             z(:,:,0)*g(0)
      cwsq_poly = cwsq_poly + z(:,:,1)*g(1)
      cwsq_poly = cwsq_poly + z(:,:,2)*g(2)

    case (1)                     ! a 2-nd order polynomial

      cwsq_poly =             z(:,:,0)*g(0)
      cwsq_poly = cwsq_poly + z(:,:,1)*g(1)

  end select
end function cwsq_poly


!    Name : cwtr_poly
! Purpose : This function computes a matrix polynomial that consists
!         : only of even order terms.
!   Input : - 2*n is the order of the polynomial.
!         : - g(i) stores the coefficient of the 2*i-th power term.
!         : - z(*,i) is the (2*i)-th power of a matrix.
!  Output : The value of the polynomial.
!    Note : This routine treats a complex, double precision, triangular matrix.

pure function cwtr_poly(n,g,z)
  integer,          intent(in) :: n
  real   (kind=wp), intent(in) :: g(0:)
  complex(kind=wp), intent(in) :: z(:,:,0:)
  complex(kind=wp)             :: cwtr_poly(size(z,1),size(z,2))
  complex(kind=wp)             :: temporary(size(z,1),size(z,2))

  integer :: j, pmax

  pmax = 0
  cwtr_poly = cmplx(0.0_wp, 0.0_wp, wp)

  ! an even polynomial of 2n-th order

  select case (n)                ! z(*,i) is (2*i)-th power of the matrix

    case (8)                     ! a 16-th order polynomial

      cwtr_poly =             z(:,:,0)*g(0)
      cwtr_poly = cwtr_poly + z(:,:,1)*g(1)
      cwtr_poly = cwtr_poly + z(:,:,2)*g(2)
      cwtr_poly = cwtr_poly + z(:,:,3)*g(3)

      temporary =             z(:,:,0)*g(4)
      temporary = temporary + z(:,:,1)*g(5)
      temporary = temporary + z(:,:,2)*g(6)
      temporary = temporary + z(:,:,3)*g(7)
      temporary = temporary + z(:,:,4)*g(8)
      temporary = temporary .trtimes. z(:,:,4)

      cwtr_poly = cwtr_poly + temporary
      pmax = 4

    case (6)                     ! a 12-th order polynomial

      cwtr_poly =             z(:,:,0)*g(0)
      cwtr_poly = cwtr_poly + z(:,:,1)*g(1)
      cwtr_poly = cwtr_poly + z(:,:,2)*g(2)

      temporary =             z(:,:,0)*g(3)
      temporary = temporary + z(:,:,1)*g(4)
      temporary = temporary + z(:,:,2)*g(5)
      temporary = temporary + z(:,:,3)*g(6)
      temporary = temporary .trtimes. z(:,:,3)

      cwtr_poly = cwtr_poly + temporary
      pmax = 3

    case (5)                     ! a 10-th order polynomial

      cwtr_poly =             z(:,:,0)*g(0)
      cwtr_poly = cwtr_poly + z(:,:,1)*g(1)
      cwtr_poly = cwtr_poly + z(:,:,2)*g(2)

      temporary =             z(:,:,0)*g(3)
      temporary = temporary + z(:,:,1)*g(4)
      temporary = temporary + z(:,:,2)*g(5)
      temporary = temporary .trtimes. z(:,:,3)

      cwtr_poly = cwtr_poly + temporary
      pmax = 3

    case (4)                     ! an 8-th order polynomial

      cwtr_poly =             z(:,:,0)*g(0)
      cwtr_poly = cwtr_poly + z(:,:,1)*g(1)

      temporary =             z(:,:,0)*g(2)
      temporary = temporary + z(:,:,1)*g(3)
      temporary = temporary + z(:,:,2)*g(4)
      temporary = temporary .trtimes. z(:,:,2)

      cwtr_poly = cwtr_poly + temporary
      pmax = 2

    case (3)                     ! a 6-th order polynomial

      cwtr_poly =             z(:,:,0)*g(0)
      cwtr_poly = cwtr_poly + z(:,:,1)*g(1)
      cwtr_poly = cwtr_poly + z(:,:,2)*g(2)
      cwtr_poly = cwtr_poly + z(:,:,3)*g(3)
      pmax = 3

    case (2)                     ! a 4-th order polynomial

      cwtr_poly =             z(:,:,0)*g(0)
      cwtr_poly = cwtr_poly + z(:,:,1)*g(1)
      cwtr_poly = cwtr_poly + z(:,:,2)*g(2)
      pmax = 2

    case (1)                     ! a 2-nd order polynomial

      cwtr_poly =             z(:,:,0)*g(0)
      cwtr_poly = cwtr_poly + z(:,:,1)*g(1)
      pmax = 1

  end select

  ! for the block structure to be invariant
  do j=1,size(z,1)-1
    if ( sum(abs(z(j+1,j,1:pmax))) < tiny(1.0_wp) ) then
      cwtr_poly(j+1,j)= cmplx(0.0_wp, 0.0_wp, wp)
    end if
  end do
end function cwtr_poly

#ifdef __USE_TPREC

!    Name : rtdg_poly
! Purpose : This function computes a matrix polynomial that consists
!         : only of even order terms.
!   Input : - 2*n is the order of the polynomial.
!         : - g(i) stores the coefficient of the 2*i-th power term.
!         : - z(*,i) is the (2*i)-th power of a matrix.
!  Output : The value of the polynomial.
!    Note : This routine treats a real, triple precision, diagonal matrix.

pure function rtdg_poly(n,g,z)
  integer,          intent(in) :: n
  real   (kind=tp), intent(in) :: g(0:)
  real   (kind=tp), intent(in) :: z(:,0:)
  real   (kind=tp)             :: rtdg_poly(size(z,1))
  real   (kind=tp)             :: temporary(size(z,1))

  rtdg_poly = 0.0_tp

  ! an even polynomial of 2n-th order

  select case (n)                ! z(*,i) is (2*i)-th power of the matrix

    case (8)                     ! a 16-th order polynomial

      rtdg_poly =             z(:,0)*g(0)
      rtdg_poly = rtdg_poly + z(:,1)*g(1)
      rtdg_poly = rtdg_poly + z(:,2)*g(2)
      rtdg_poly = rtdg_poly + z(:,3)*g(3)

      temporary =             z(:,0)*g(4)
      temporary = temporary + z(:,1)*g(5)
      temporary = temporary + z(:,2)*g(6)
      temporary = temporary + z(:,3)*g(7)
      temporary = temporary + z(:,4)*g(8)
      temporary = temporary .dgtimes. z(:,4)

      rtdg_poly = rtdg_poly + temporary

    case (6)                     ! a 12-th order polynomial

      rtdg_poly =             z(:,0)*g(0)
      rtdg_poly = rtdg_poly + z(:,1)*g(1)
      rtdg_poly = rtdg_poly + z(:,2)*g(2)

      temporary =             z(:,0)*g(3)
      temporary = temporary + z(:,1)*g(4)
      temporary = temporary + z(:,2)*g(5)
      temporary = temporary + z(:,3)*g(6)
      temporary = temporary .dgtimes. z(:,3)

      rtdg_poly = rtdg_poly + temporary

    case (5)                     ! a 10-th order polynomial

      rtdg_poly =             z(:,0)*g(0)
      rtdg_poly = rtdg_poly + z(:,1)*g(1)
      rtdg_poly = rtdg_poly + z(:,2)*g(2)

      temporary =             z(:,0)*g(3)
      temporary = temporary + z(:,1)*g(4)
      temporary = temporary + z(:,2)*g(5)
      temporary = temporary .dgtimes. z(:,3)

      rtdg_poly = rtdg_poly + temporary

    case (4)                     ! an 8-th order polynomial

      rtdg_poly =             z(:,0)*g(0)
      rtdg_poly = rtdg_poly + z(:,1)*g(1)

      temporary =             z(:,0)*g(2)
      temporary = temporary + z(:,1)*g(3)
      temporary = temporary + z(:,2)*g(4)
      temporary = temporary .dgtimes. z(:,2)

      rtdg_poly = rtdg_poly + temporary

    case (3)                     ! a 6-th order polynomial

      rtdg_poly =             z(:,0)*g(0)
      rtdg_poly = rtdg_poly + z(:,1)*g(1)
      rtdg_poly = rtdg_poly + z(:,2)*g(2)
      rtdg_poly = rtdg_poly + z(:,3)*g(3)

    case (2)                     ! a 4-th order polynomial

      rtdg_poly =             z(:,0)*g(0)
      rtdg_poly = rtdg_poly + z(:,1)*g(1)
      rtdg_poly = rtdg_poly + z(:,2)*g(2)

    case (1)                     ! a 2-nd order polynomial

      rtdg_poly =             z(:,0)*g(0)
      rtdg_poly = rtdg_poly + z(:,1)*g(1)

  end select
end function rtdg_poly


!    Name : rtsq_poly
! Purpose : This function computes a matrix polynomial that consists
!         : only of even order terms.
!   Input : - 2*n is the order of the polynomial.
!         : - g(i) stores the coefficient of the 2*i-th power term.
!         : - z(*,i) is the (2*i)-th power of a matrix.
!  Output : The value of the polynomial.
!    Note : This routine treats a real, triple precision, square matrix.

pure function rtsq_poly(n,g,z)
  integer,          intent(in) :: n
  real   (kind=tp), intent(in) :: g(0:)
  real   (kind=tp), intent(in) :: z(:,:,0:)
  real   (kind=tp)             :: rtsq_poly(size(z,1),size(z,2))
  real   (kind=tp)             :: temporary(size(z,1),size(z,2))

  rtsq_poly = 0.0_tp

  ! an even polynomial of 2n-th order

  select case (n)                ! z(*,i) is (2*i)-th power of the matrix

    case (8)                     ! a 16-th order polynomial

      rtsq_poly =             z(:,:,0)*g(0)
      rtsq_poly = rtsq_poly + z(:,:,1)*g(1)
      rtsq_poly = rtsq_poly + z(:,:,2)*g(2)
      rtsq_poly = rtsq_poly + z(:,:,3)*g(3)

      temporary =             z(:,:,0)*g(4)
      temporary = temporary + z(:,:,1)*g(5)
      temporary = temporary + z(:,:,2)*g(6)
      temporary = temporary + z(:,:,3)*g(7)
      temporary = temporary + z(:,:,4)*g(8)
      temporary = temporary .sqtimes. z(:,:,4)

      rtsq_poly = rtsq_poly + temporary

    case (6)                     ! a 12-th order polynomial

      rtsq_poly =             z(:,:,0)*g(0)
      rtsq_poly = rtsq_poly + z(:,:,1)*g(1)
      rtsq_poly = rtsq_poly + z(:,:,2)*g(2)

      temporary =             z(:,:,0)*g(3)
      temporary = temporary + z(:,:,1)*g(4)
      temporary = temporary + z(:,:,2)*g(5)
      temporary = temporary + z(:,:,3)*g(6)
      temporary = temporary .sqtimes. z(:,:,3)

      rtsq_poly = rtsq_poly + temporary

    case (5)                     ! a 10-th order polynomial

      rtsq_poly =             z(:,:,0)*g(0)
      rtsq_poly = rtsq_poly + z(:,:,1)*g(1)
      rtsq_poly = rtsq_poly + z(:,:,2)*g(2)

      temporary =             z(:,:,0)*g(3)
      temporary = temporary + z(:,:,1)*g(4)
      temporary = temporary + z(:,:,2)*g(5)
      temporary = temporary .sqtimes. z(:,:,3)

      rtsq_poly = rtsq_poly + temporary

    case (4)                     ! an 8-th order polynomial

      rtsq_poly =             z(:,:,0)*g(0)
      rtsq_poly = rtsq_poly + z(:,:,1)*g(1)

      temporary =             z(:,:,0)*g(2)
      temporary = temporary + z(:,:,1)*g(3)
      temporary = temporary + z(:,:,2)*g(4)
      temporary = temporary .sqtimes. z(:,:,2)

      rtsq_poly = rtsq_poly + temporary

    case (3)                     ! a 6-th order polynomial

      rtsq_poly =             z(:,:,0)*g(0)
      rtsq_poly = rtsq_poly + z(:,:,1)*g(1)
      rtsq_poly = rtsq_poly + z(:,:,2)*g(2)
      rtsq_poly = rtsq_poly + z(:,:,3)*g(3)

    case (2)                     ! a 4-th order polynomial

      rtsq_poly =             z(:,:,0)*g(0)
      rtsq_poly = rtsq_poly + z(:,:,1)*g(1)
      rtsq_poly = rtsq_poly + z(:,:,2)*g(2)

    case (1)                     ! a 2-nd order polynomial

      rtsq_poly =             z(:,:,0)*g(0)
      rtsq_poly = rtsq_poly + z(:,:,1)*g(1)

  end select
end function rtsq_poly


!    Name : rttr_poly
! Purpose : This function computes a matrix polynomial that consists
!         : only of even order terms.
!   Input : - 2*n is the order of the polynomial.
!         : - g(i) stores the coefficient of the 2*i-th power term.
!         : - z(*,i) is the (2*i)-th power of a matrix.
!  Output : The value of the polynomial.
!    Note : This routine treats a real, triple precision, triangular matrix.

pure function rttr_poly(n,g,z)
  integer,          intent(in) :: n
  real   (kind=tp), intent(in) :: g(0:)
  real   (kind=tp), intent(in) :: z(:,:,0:)
  real   (kind=tp)             :: rttr_poly(size(z,1),size(z,2))
  real   (kind=tp)             :: temporary(size(z,1),size(z,2))

  integer :: j, pmax

  pmax = 0
  rttr_poly = 0.0_tp

  ! an even polynomial of 2n-th order

  select case (n)                ! z(*,i) is (2*i)-th power of the matrix

    case (8)                     ! a 16-th order polynomial

      rttr_poly =             z(:,:,0)*g(0)
      rttr_poly = rttr_poly + z(:,:,1)*g(1)
      rttr_poly = rttr_poly + z(:,:,2)*g(2)
      rttr_poly = rttr_poly + z(:,:,3)*g(3)

      temporary =             z(:,:,0)*g(4)
      temporary = temporary + z(:,:,1)*g(5)
      temporary = temporary + z(:,:,2)*g(6)
      temporary = temporary + z(:,:,3)*g(7)
      temporary = temporary + z(:,:,4)*g(8)
      temporary = temporary .trtimes. z(:,:,4)

      rttr_poly = rttr_poly + temporary
      pmax = 4

    case (6)                     ! a 12-th order polynomial

      rttr_poly =             z(:,:,0)*g(0)
      rttr_poly = rttr_poly + z(:,:,1)*g(1)
      rttr_poly = rttr_poly + z(:,:,2)*g(2)

      temporary =             z(:,:,0)*g(3)
      temporary = temporary + z(:,:,1)*g(4)
      temporary = temporary + z(:,:,2)*g(5)
      temporary = temporary + z(:,:,3)*g(6)
      temporary = temporary .trtimes. z(:,:,3)

      rttr_poly = rttr_poly + temporary
      pmax = 3

    case (5)                     ! a 10-th order polynomial

      rttr_poly =             z(:,:,0)*g(0)
      rttr_poly = rttr_poly + z(:,:,1)*g(1)
      rttr_poly = rttr_poly + z(:,:,2)*g(2)

      temporary =             z(:,:,0)*g(3)
      temporary = temporary + z(:,:,1)*g(4)
      temporary = temporary + z(:,:,2)*g(5)
      temporary = temporary .trtimes. z(:,:,3)

      rttr_poly = rttr_poly + temporary
      pmax = 3

    case (4)                     ! an 8-th order polynomial

      rttr_poly =             z(:,:,0)*g(0)
      rttr_poly = rttr_poly + z(:,:,1)*g(1)

      temporary =             z(:,:,0)*g(2)
      temporary = temporary + z(:,:,1)*g(3)
      temporary = temporary + z(:,:,2)*g(4)
      temporary = temporary .trtimes. z(:,:,2)

      rttr_poly = rttr_poly + temporary
      pmax = 2

    case (3)                     ! a 6-th order polynomial

      rttr_poly =             z(:,:,0)*g(0)
      rttr_poly = rttr_poly + z(:,:,1)*g(1)
      rttr_poly = rttr_poly + z(:,:,2)*g(2)
      rttr_poly = rttr_poly + z(:,:,3)*g(3)
      pmax = 3

    case (2)                     ! a 4-th order polynomial

      rttr_poly =             z(:,:,0)*g(0)
      rttr_poly = rttr_poly + z(:,:,1)*g(1)
      rttr_poly = rttr_poly + z(:,:,2)*g(2)
      pmax = 2

    case (1)                     ! a 2-nd order polynomial

      rttr_poly =             z(:,:,0)*g(0)
      rttr_poly = rttr_poly + z(:,:,1)*g(1)
      pmax = 1

  end select

  ! for the block structure to be invariant
  do j=1,size(z,1)-1
    if ( sum(abs(z(j+1,j,1:pmax))) < tiny(1.0_tp) ) then
      rttr_poly(j+1,j)= 0.0_tp
    end if
  end do
end function rttr_poly


!    Name : ctdg_poly
! Purpose : This function computes a matrix polynomial that consists
!         : only of even order terms.
!   Input : - 2*n is the order of the polynomial.
!         : - g(i) stores the coefficient of the 2*i-th power term.
!         : - z(*,i) is the (2*i)-th power of a matrix.
!  Output : The value of the polynomial.
!    Note : This routine treats a complex, triple precision, diagonal matrix.

pure function ctdg_poly(n,g,z)
  integer,          intent(in) :: n
  real   (kind=tp), intent(in) :: g(0:)
  complex(kind=tp), intent(in) :: z(:,0:)
  complex(kind=tp)             :: ctdg_poly(size(z,1))
  complex(kind=tp)             :: temporary(size(z,1))

  ctdg_poly = cmplx(0.0_tp, 0.0_tp, tp)

  ! an even polynomial of 2n-th order

  select case (n)                ! z(*,i) is (2*i)-th power of the matrix

    case (8)                     ! a 16-th order polynomial

      ctdg_poly =             z(:,0)*g(0)
      ctdg_poly = ctdg_poly + z(:,1)*g(1)
      ctdg_poly = ctdg_poly + z(:,2)*g(2)
      ctdg_poly = ctdg_poly + z(:,3)*g(3)

      temporary =             z(:,0)*g(4)
      temporary = temporary + z(:,1)*g(5)
      temporary = temporary + z(:,2)*g(6)
      temporary = temporary + z(:,3)*g(7)
      temporary = temporary + z(:,4)*g(8)
      temporary = temporary .dgtimes. z(:,4)

      ctdg_poly = ctdg_poly + temporary

    case (6)                     ! a 12-th order polynomial

      ctdg_poly =             z(:,0)*g(0)
      ctdg_poly = ctdg_poly + z(:,1)*g(1)
      ctdg_poly = ctdg_poly + z(:,2)*g(2)

      temporary =             z(:,0)*g(3)
      temporary = temporary + z(:,1)*g(4)
      temporary = temporary + z(:,2)*g(5)
      temporary = temporary + z(:,3)*g(6)
      temporary = temporary .dgtimes. z(:,3)

      ctdg_poly = ctdg_poly + temporary

    case (5)                     ! a 10-th order polynomial

      ctdg_poly =             z(:,0)*g(0)
      ctdg_poly = ctdg_poly + z(:,1)*g(1)
      ctdg_poly = ctdg_poly + z(:,2)*g(2)

      temporary =             z(:,0)*g(3)
      temporary = temporary + z(:,1)*g(4)
      temporary = temporary + z(:,2)*g(5)
      temporary = temporary .dgtimes. z(:,3)

      ctdg_poly = ctdg_poly + temporary

    case (4)                     ! an 8-th order polynomial

      ctdg_poly =             z(:,0)*g(0)
      ctdg_poly = ctdg_poly + z(:,1)*g(1)

      temporary =             z(:,0)*g(2)
      temporary = temporary + z(:,1)*g(3)
      temporary = temporary + z(:,2)*g(4)
      temporary = temporary .dgtimes. z(:,2)

      ctdg_poly = ctdg_poly + temporary

    case (3)                     ! a 6-th order polynomial

      ctdg_poly =             z(:,0)*g(0)
      ctdg_poly = ctdg_poly + z(:,1)*g(1)
      ctdg_poly = ctdg_poly + z(:,2)*g(2)
      ctdg_poly = ctdg_poly + z(:,3)*g(3)

    case (2)                     ! a 4-th order polynomial

      ctdg_poly =             z(:,0)*g(0)
      ctdg_poly = ctdg_poly + z(:,1)*g(1)
      ctdg_poly = ctdg_poly + z(:,2)*g(2)

    case (1)                     ! a 2-nd order polynomial

      ctdg_poly =             z(:,0)*g(0)
      ctdg_poly = ctdg_poly + z(:,1)*g(1)

  end select
end function ctdg_poly


!    Name : ctsq_poly
! Purpose : This function computes a matrix polynomial that consists
!         : only of even order terms.
!   Input : - 2*n is the order of the polynomial.
!         : - g(i) stores the coefficient of the 2*i-th power term.
!         : - z(*,i) is the (2*i)-th power of a matrix.
!  Output : The value of the polynomial.
!    Note : This routine treats a complex, triple precision, square matrix.

pure function ctsq_poly(n,g,z)
  integer,          intent(in) :: n
  real   (kind=tp), intent(in) :: g(0:)
  complex(kind=tp), intent(in) :: z(:,:,0:)
  complex(kind=tp)             :: ctsq_poly(size(z,1),size(z,2))
  complex(kind=tp)             :: temporary(size(z,1),size(z,2))

  ctsq_poly = cmplx(0.0_tp, 0.0_tp, tp)

  ! an even polynomial of 2n-th order

  select case (n)                ! z(*,i) is (2*i)-th power of the matrix

    case (8)                     ! a 16-th order polynomial

      ctsq_poly =             z(:,:,0)*g(0)
      ctsq_poly = ctsq_poly + z(:,:,1)*g(1)
      ctsq_poly = ctsq_poly + z(:,:,2)*g(2)
      ctsq_poly = ctsq_poly + z(:,:,3)*g(3)

      temporary =             z(:,:,0)*g(4)
      temporary = temporary + z(:,:,1)*g(5)
      temporary = temporary + z(:,:,2)*g(6)
      temporary = temporary + z(:,:,3)*g(7)
      temporary = temporary + z(:,:,4)*g(8)
      temporary = temporary .sqtimes. z(:,:,4)

      ctsq_poly = ctsq_poly + temporary

    case (6)                     ! a 12-th order polynomial

      ctsq_poly =             z(:,:,0)*g(0)
      ctsq_poly = ctsq_poly + z(:,:,1)*g(1)
      ctsq_poly = ctsq_poly + z(:,:,2)*g(2)

      temporary =             z(:,:,0)*g(3)
      temporary = temporary + z(:,:,1)*g(4)
      temporary = temporary + z(:,:,2)*g(5)
      temporary = temporary + z(:,:,3)*g(6)
      temporary = temporary .sqtimes. z(:,:,3)

      ctsq_poly = ctsq_poly + temporary

    case (5)                     ! a 10-th order polynomial

      ctsq_poly =             z(:,:,0)*g(0)
      ctsq_poly = ctsq_poly + z(:,:,1)*g(1)
      ctsq_poly = ctsq_poly + z(:,:,2)*g(2)

      temporary =             z(:,:,0)*g(3)
      temporary = temporary + z(:,:,1)*g(4)
      temporary = temporary + z(:,:,2)*g(5)
      temporary = temporary .sqtimes. z(:,:,3)

      ctsq_poly = ctsq_poly + temporary

    case (4)                     ! an 8-th order polynomial

      ctsq_poly =             z(:,:,0)*g(0)
      ctsq_poly = ctsq_poly + z(:,:,1)*g(1)

      temporary =             z(:,:,0)*g(2)
      temporary = temporary + z(:,:,1)*g(3)
      temporary = temporary + z(:,:,2)*g(4)
      temporary = temporary .sqtimes. z(:,:,2)

      ctsq_poly = ctsq_poly + temporary

    case (3)                     ! a 6-th order polynomial

      ctsq_poly =             z(:,:,0)*g(0)
      ctsq_poly = ctsq_poly + z(:,:,1)*g(1)
      ctsq_poly = ctsq_poly + z(:,:,2)*g(2)
      ctsq_poly = ctsq_poly + z(:,:,3)*g(3)

    case (2)                     ! a 4-th order polynomial

      ctsq_poly =             z(:,:,0)*g(0)
      ctsq_poly = ctsq_poly + z(:,:,1)*g(1)
      ctsq_poly = ctsq_poly + z(:,:,2)*g(2)

    case (1)                     ! a 2-nd order polynomial

      ctsq_poly =             z(:,:,0)*g(0)
      ctsq_poly = ctsq_poly + z(:,:,1)*g(1)

  end select
end function ctsq_poly


!    Name : cttr_poly
! Purpose : This function computes a matrix polynomial that consists
!         : only of even order terms.
!   Input : - 2*n is the order of the polynomial.
!         : - g(i) stores the coefficient of the 2*i-th power term.
!         : - z(*,i) is the (2*i)-th power of a matrix.
!  Output : The value of the polynomial.
!    Note : This routine treats a complex, triple precision, triangular matrix.

pure function cttr_poly(n,g,z)
  integer,          intent(in) :: n
  real   (kind=tp), intent(in) :: g(0:)
  complex(kind=tp), intent(in) :: z(:,:,0:)
  complex(kind=tp)             :: cttr_poly(size(z,1),size(z,2))
  complex(kind=tp)             :: temporary(size(z,1),size(z,2))

  integer :: j, pmax

  pmax = 0
  cttr_poly = cmplx(0.0_tp, 0.0_tp, tp)

  ! an even polynomial of 2n-th order

  select case (n)                ! z(*,i) is (2*i)-th power of the matrix

    case (8)                     ! a 16-th order polynomial

      cttr_poly =             z(:,:,0)*g(0)
      cttr_poly = cttr_poly + z(:,:,1)*g(1)
      cttr_poly = cttr_poly + z(:,:,2)*g(2)
      cttr_poly = cttr_poly + z(:,:,3)*g(3)

      temporary =             z(:,:,0)*g(4)
      temporary = temporary + z(:,:,1)*g(5)
      temporary = temporary + z(:,:,2)*g(6)
      temporary = temporary + z(:,:,3)*g(7)
      temporary = temporary + z(:,:,4)*g(8)
      temporary = temporary .trtimes. z(:,:,4)

      cttr_poly = cttr_poly + temporary
      pmax = 4

    case (6)                     ! a 12-th order polynomial

      cttr_poly =             z(:,:,0)*g(0)
      cttr_poly = cttr_poly + z(:,:,1)*g(1)
      cttr_poly = cttr_poly + z(:,:,2)*g(2)

      temporary =             z(:,:,0)*g(3)
      temporary = temporary + z(:,:,1)*g(4)
      temporary = temporary + z(:,:,2)*g(5)
      temporary = temporary + z(:,:,3)*g(6)
      temporary = temporary .trtimes. z(:,:,3)

      cttr_poly = cttr_poly + temporary
      pmax = 3

    case (5)                     ! a 10-th order polynomial

      cttr_poly =             z(:,:,0)*g(0)
      cttr_poly = cttr_poly + z(:,:,1)*g(1)
      cttr_poly = cttr_poly + z(:,:,2)*g(2)

      temporary =             z(:,:,0)*g(3)
      temporary = temporary + z(:,:,1)*g(4)
      temporary = temporary + z(:,:,2)*g(5)
      temporary = temporary .trtimes. z(:,:,3)

      cttr_poly = cttr_poly + temporary
      pmax = 3

    case (4)                     ! an 8-th order polynomial

      cttr_poly =             z(:,:,0)*g(0)
      cttr_poly = cttr_poly + z(:,:,1)*g(1)

      temporary =             z(:,:,0)*g(2)
      temporary = temporary + z(:,:,1)*g(3)
      temporary = temporary + z(:,:,2)*g(4)
      temporary = temporary .trtimes. z(:,:,2)

      cttr_poly = cttr_poly + temporary
      pmax = 2

    case (3)                     ! a 6-th order polynomial

      cttr_poly =             z(:,:,0)*g(0)
      cttr_poly = cttr_poly + z(:,:,1)*g(1)
      cttr_poly = cttr_poly + z(:,:,2)*g(2)
      cttr_poly = cttr_poly + z(:,:,3)*g(3)
      pmax = 3

    case (2)                     ! a 4-th order polynomial

      cttr_poly =             z(:,:,0)*g(0)
      cttr_poly = cttr_poly + z(:,:,1)*g(1)
      cttr_poly = cttr_poly + z(:,:,2)*g(2)
      pmax = 2

    case (1)                     ! a 2-nd order polynomial

      cttr_poly =             z(:,:,0)*g(0)
      cttr_poly = cttr_poly + z(:,:,1)*g(1)
      pmax = 1

  end select

  ! for the block structure to be invariant
  do j=1,size(z,1)-1
    if ( sum(abs(z(j+1,j,1:pmax))) < tiny(1.0_tp) ) then
      cttr_poly(j+1,j)= cmplx(0.0_tp, 0.0_tp, tp)
    end if
  end do
end function cttr_poly

#endif
#ifdef __USE_QPREC

!    Name : rqdg_poly
! Purpose : This function computes a matrix polynomial that consists
!         : only of even order terms.
!   Input : - 2*n is the order of the polynomial.
!         : - g(i) stores the coefficient of the 2*i-th power term.
!         : - z(*,i) is the (2*i)-th power of a matrix.
!  Output : The value of the polynomial.
!    Note : This routine treats a real, quadruple precision, diagonal matrix.

pure function rqdg_poly(n,g,z)
  integer,          intent(in) :: n
  real   (kind=qp), intent(in) :: g(0:)
  real   (kind=qp), intent(in) :: z(:,0:)
  real   (kind=qp)             :: rqdg_poly(size(z,1))
  real   (kind=qp)             :: temporary(size(z,1))

  rqdg_poly = 0.0_qp

  ! an even polynomial of 2n-th order

  select case (n)                ! z(*,i) is (2*i)-th power of the matrix

    case (8)                     ! a 16-th order polynomial

      rqdg_poly =             z(:,0)*g(0)
      rqdg_poly = rqdg_poly + z(:,1)*g(1)
      rqdg_poly = rqdg_poly + z(:,2)*g(2)
      rqdg_poly = rqdg_poly + z(:,3)*g(3)

      temporary =             z(:,0)*g(4)
      temporary = temporary + z(:,1)*g(5)
      temporary = temporary + z(:,2)*g(6)
      temporary = temporary + z(:,3)*g(7)
      temporary = temporary + z(:,4)*g(8)
      temporary = temporary .dgtimes. z(:,4)

      rqdg_poly = rqdg_poly + temporary

    case (6)                     ! a 12-th order polynomial

      rqdg_poly =             z(:,0)*g(0)
      rqdg_poly = rqdg_poly + z(:,1)*g(1)
      rqdg_poly = rqdg_poly + z(:,2)*g(2)

      temporary =             z(:,0)*g(3)
      temporary = temporary + z(:,1)*g(4)
      temporary = temporary + z(:,2)*g(5)
      temporary = temporary + z(:,3)*g(6)
      temporary = temporary .dgtimes. z(:,3)

      rqdg_poly = rqdg_poly + temporary

    case (5)                     ! a 10-th order polynomial

      rqdg_poly =             z(:,0)*g(0)
      rqdg_poly = rqdg_poly + z(:,1)*g(1)
      rqdg_poly = rqdg_poly + z(:,2)*g(2)

      temporary =             z(:,0)*g(3)
      temporary = temporary + z(:,1)*g(4)
      temporary = temporary + z(:,2)*g(5)
      temporary = temporary .dgtimes. z(:,3)

      rqdg_poly = rqdg_poly + temporary

    case (4)                     ! an 8-th order polynomial

      rqdg_poly =             z(:,0)*g(0)
      rqdg_poly = rqdg_poly + z(:,1)*g(1)

      temporary =             z(:,0)*g(2)
      temporary = temporary + z(:,1)*g(3)
      temporary = temporary + z(:,2)*g(4)
      temporary = temporary .dgtimes. z(:,2)

      rqdg_poly = rqdg_poly + temporary

    case (3)                     ! a 6-th order polynomial

      rqdg_poly =             z(:,0)*g(0)
      rqdg_poly = rqdg_poly + z(:,1)*g(1)
      rqdg_poly = rqdg_poly + z(:,2)*g(2)
      rqdg_poly = rqdg_poly + z(:,3)*g(3)

    case (2)                     ! a 4-th order polynomial

      rqdg_poly =             z(:,0)*g(0)
      rqdg_poly = rqdg_poly + z(:,1)*g(1)
      rqdg_poly = rqdg_poly + z(:,2)*g(2)

    case (1)                     ! a 2-nd order polynomial

      rqdg_poly =             z(:,0)*g(0)
      rqdg_poly = rqdg_poly + z(:,1)*g(1)

  end select
end function rqdg_poly


!    Name : rqsq_poly
! Purpose : This function computes a matrix polynomial that consists
!         : only of even order terms.
!   Input : - 2*n is the order of the polynomial.
!         : - g(i) stores the coefficient of the 2*i-th power term.
!         : - z(*,i) is the (2*i)-th power of a matrix.
!  Output : The value of the polynomial.
!    Note : This routine treats a real, quadruple precision, square matrix.

pure function rqsq_poly(n,g,z)
  integer,          intent(in) :: n
  real   (kind=qp), intent(in) :: g(0:)
  real   (kind=qp), intent(in) :: z(:,:,0:)
  real   (kind=qp)             :: rqsq_poly(size(z,1),size(z,2))
  real   (kind=qp)             :: temporary(size(z,1),size(z,2))

  rqsq_poly = 0.0_qp

  ! an even polynomial of 2n-th order

  select case (n)                ! z(*,i) is (2*i)-th power of the matrix

    case (8)                     ! a 16-th order polynomial

      rqsq_poly =             z(:,:,0)*g(0)
      rqsq_poly = rqsq_poly + z(:,:,1)*g(1)
      rqsq_poly = rqsq_poly + z(:,:,2)*g(2)
      rqsq_poly = rqsq_poly + z(:,:,3)*g(3)

      temporary =             z(:,:,0)*g(4)
      temporary = temporary + z(:,:,1)*g(5)
      temporary = temporary + z(:,:,2)*g(6)
      temporary = temporary + z(:,:,3)*g(7)
      temporary = temporary + z(:,:,4)*g(8)
      temporary = temporary .sqtimes. z(:,:,4)

      rqsq_poly = rqsq_poly + temporary

    case (6)                     ! a 12-th order polynomial

      rqsq_poly =             z(:,:,0)*g(0)
      rqsq_poly = rqsq_poly + z(:,:,1)*g(1)
      rqsq_poly = rqsq_poly + z(:,:,2)*g(2)

      temporary =             z(:,:,0)*g(3)
      temporary = temporary + z(:,:,1)*g(4)
      temporary = temporary + z(:,:,2)*g(5)
      temporary = temporary + z(:,:,3)*g(6)
      temporary = temporary .sqtimes. z(:,:,3)

      rqsq_poly = rqsq_poly + temporary

    case (5)                     ! a 10-th order polynomial

      rqsq_poly =             z(:,:,0)*g(0)
      rqsq_poly = rqsq_poly + z(:,:,1)*g(1)
      rqsq_poly = rqsq_poly + z(:,:,2)*g(2)

      temporary =             z(:,:,0)*g(3)
      temporary = temporary + z(:,:,1)*g(4)
      temporary = temporary + z(:,:,2)*g(5)
      temporary = temporary .sqtimes. z(:,:,3)

      rqsq_poly = rqsq_poly + temporary

    case (4)                     ! an 8-th order polynomial

      rqsq_poly =             z(:,:,0)*g(0)
      rqsq_poly = rqsq_poly + z(:,:,1)*g(1)

      temporary =             z(:,:,0)*g(2)
      temporary = temporary + z(:,:,1)*g(3)
      temporary = temporary + z(:,:,2)*g(4)
      temporary = temporary .sqtimes. z(:,:,2)

      rqsq_poly = rqsq_poly + temporary

    case (3)                     ! a 6-th order polynomial

      rqsq_poly =             z(:,:,0)*g(0)
      rqsq_poly = rqsq_poly + z(:,:,1)*g(1)
      rqsq_poly = rqsq_poly + z(:,:,2)*g(2)
      rqsq_poly = rqsq_poly + z(:,:,3)*g(3)

    case (2)                     ! a 4-th order polynomial

      rqsq_poly =             z(:,:,0)*g(0)
      rqsq_poly = rqsq_poly + z(:,:,1)*g(1)
      rqsq_poly = rqsq_poly + z(:,:,2)*g(2)

    case (1)                     ! a 2-nd order polynomial

      rqsq_poly =             z(:,:,0)*g(0)
      rqsq_poly = rqsq_poly + z(:,:,1)*g(1)

  end select
end function rqsq_poly


!    Name : rqtr_poly
! Purpose : This function computes a matrix polynomial that consists
!         : only of even order terms.
!   Input : - 2*n is the order of the polynomial.
!         : - g(i) stores the coefficient of the 2*i-th power term.
!         : - z(*,i) is the (2*i)-th power of a matrix.
!  Output : The value of the polynomial.
!    Note : This routine treats a real, quadruple precision, triangular matrix.

pure function rqtr_poly(n,g,z)
  integer,          intent(in) :: n
  real   (kind=qp), intent(in) :: g(0:)
  real   (kind=qp), intent(in) :: z(:,:,0:)
  real   (kind=qp)             :: rqtr_poly(size(z,1),size(z,2))
  real   (kind=qp)             :: temporary(size(z,1),size(z,2))

  integer :: j, pmax

  pmax = 0
  rqtr_poly = 0.0_qp

  ! an even polynomial of 2n-th order

  select case (n)                ! z(*,i) is (2*i)-th power of the matrix

    case (8)                     ! a 16-th order polynomial

      rqtr_poly =             z(:,:,0)*g(0)
      rqtr_poly = rqtr_poly + z(:,:,1)*g(1)
      rqtr_poly = rqtr_poly + z(:,:,2)*g(2)
      rqtr_poly = rqtr_poly + z(:,:,3)*g(3)

      temporary =             z(:,:,0)*g(4)
      temporary = temporary + z(:,:,1)*g(5)
      temporary = temporary + z(:,:,2)*g(6)
      temporary = temporary + z(:,:,3)*g(7)
      temporary = temporary + z(:,:,4)*g(8)
      temporary = temporary .trtimes. z(:,:,4)

      rqtr_poly = rqtr_poly + temporary
      pmax = 4

    case (6)                     ! a 12-th order polynomial

      rqtr_poly =             z(:,:,0)*g(0)
      rqtr_poly = rqtr_poly + z(:,:,1)*g(1)
      rqtr_poly = rqtr_poly + z(:,:,2)*g(2)

      temporary =             z(:,:,0)*g(3)
      temporary = temporary + z(:,:,1)*g(4)
      temporary = temporary + z(:,:,2)*g(5)
      temporary = temporary + z(:,:,3)*g(6)
      temporary = temporary .trtimes. z(:,:,3)

      rqtr_poly = rqtr_poly + temporary
      pmax = 3

    case (5)                     ! a 10-th order polynomial

      rqtr_poly =             z(:,:,0)*g(0)
      rqtr_poly = rqtr_poly + z(:,:,1)*g(1)
      rqtr_poly = rqtr_poly + z(:,:,2)*g(2)

      temporary =             z(:,:,0)*g(3)
      temporary = temporary + z(:,:,1)*g(4)
      temporary = temporary + z(:,:,2)*g(5)
      temporary = temporary .trtimes. z(:,:,3)

      rqtr_poly = rqtr_poly + temporary
      pmax = 3

    case (4)                     ! an 8-th order polynomial

      rqtr_poly =             z(:,:,0)*g(0)
      rqtr_poly = rqtr_poly + z(:,:,1)*g(1)

      temporary =             z(:,:,0)*g(2)
      temporary = temporary + z(:,:,1)*g(3)
      temporary = temporary + z(:,:,2)*g(4)
      temporary = temporary .trtimes. z(:,:,2)

      rqtr_poly = rqtr_poly + temporary
      pmax = 2

    case (3)                     ! a 6-th order polynomial

      rqtr_poly =             z(:,:,0)*g(0)
      rqtr_poly = rqtr_poly + z(:,:,1)*g(1)
      rqtr_poly = rqtr_poly + z(:,:,2)*g(2)
      rqtr_poly = rqtr_poly + z(:,:,3)*g(3)
      pmax = 3

    case (2)                     ! a 4-th order polynomial

      rqtr_poly =             z(:,:,0)*g(0)
      rqtr_poly = rqtr_poly + z(:,:,1)*g(1)
      rqtr_poly = rqtr_poly + z(:,:,2)*g(2)
      pmax = 2

    case (1)                     ! a 2-nd order polynomial

      rqtr_poly =             z(:,:,0)*g(0)
      rqtr_poly = rqtr_poly + z(:,:,1)*g(1)
      pmax = 1

  end select

  ! for the block structure to be invariant
  do j=1,size(z,1)-1
    if ( sum(abs(z(j+1,j,1:pmax))) < tiny(1.0_qp) ) then
      rqtr_poly(j+1,j)= 0.0_qp
    end if
  end do
end function rqtr_poly


!    Name : cqdg_poly
! Purpose : This function computes a matrix polynomial that consists
!         : only of even order terms.
!   Input : - 2*n is the order of the polynomial.
!         : - g(i) stores the coefficient of the 2*i-th power term.
!         : - z(*,i) is the (2*i)-th power of a matrix.
!  Output : The value of the polynomial.
!    Note : This routine treats a complex, quadruple precision, diagonal matrix.

pure function cqdg_poly(n,g,z)
  integer,          intent(in) :: n
  real   (kind=qp), intent(in) :: g(0:)
  complex(kind=qp), intent(in) :: z(:,0:)
  complex(kind=qp)             :: cqdg_poly(size(z,1))
  complex(kind=qp)             :: temporary(size(z,1))

  cqdg_poly = cmplx(0.0_qp, 0.0_qp, qp)

  ! an even polynomial of 2n-th order

  select case (n)                ! z(*,i) is (2*i)-th power of the matrix

    case (8)                     ! a 16-th order polynomial

      cqdg_poly =             z(:,0)*g(0)
      cqdg_poly = cqdg_poly + z(:,1)*g(1)
      cqdg_poly = cqdg_poly + z(:,2)*g(2)
      cqdg_poly = cqdg_poly + z(:,3)*g(3)

      temporary =             z(:,0)*g(4)
      temporary = temporary + z(:,1)*g(5)
      temporary = temporary + z(:,2)*g(6)
      temporary = temporary + z(:,3)*g(7)
      temporary = temporary + z(:,4)*g(8)
      temporary = temporary .dgtimes. z(:,4)

      cqdg_poly = cqdg_poly + temporary

    case (6)                     ! a 12-th order polynomial

      cqdg_poly =             z(:,0)*g(0)
      cqdg_poly = cqdg_poly + z(:,1)*g(1)
      cqdg_poly = cqdg_poly + z(:,2)*g(2)

      temporary =             z(:,0)*g(3)
      temporary = temporary + z(:,1)*g(4)
      temporary = temporary + z(:,2)*g(5)
      temporary = temporary + z(:,3)*g(6)
      temporary = temporary .dgtimes. z(:,3)

      cqdg_poly = cqdg_poly + temporary

    case (5)                     ! a 10-th order polynomial

      cqdg_poly =             z(:,0)*g(0)
      cqdg_poly = cqdg_poly + z(:,1)*g(1)
      cqdg_poly = cqdg_poly + z(:,2)*g(2)

      temporary =             z(:,0)*g(3)
      temporary = temporary + z(:,1)*g(4)
      temporary = temporary + z(:,2)*g(5)
      temporary = temporary .dgtimes. z(:,3)

      cqdg_poly = cqdg_poly + temporary

    case (4)                     ! an 8-th order polynomial

      cqdg_poly =             z(:,0)*g(0)
      cqdg_poly = cqdg_poly + z(:,1)*g(1)

      temporary =             z(:,0)*g(2)
      temporary = temporary + z(:,1)*g(3)
      temporary = temporary + z(:,2)*g(4)
      temporary = temporary .dgtimes. z(:,2)

      cqdg_poly = cqdg_poly + temporary

    case (3)                     ! a 6-th order polynomial

      cqdg_poly =             z(:,0)*g(0)
      cqdg_poly = cqdg_poly + z(:,1)*g(1)
      cqdg_poly = cqdg_poly + z(:,2)*g(2)
      cqdg_poly = cqdg_poly + z(:,3)*g(3)

    case (2)                     ! a 4-th order polynomial

      cqdg_poly =             z(:,0)*g(0)
      cqdg_poly = cqdg_poly + z(:,1)*g(1)
      cqdg_poly = cqdg_poly + z(:,2)*g(2)

    case (1)                     ! a 2-nd order polynomial

      cqdg_poly =             z(:,0)*g(0)
      cqdg_poly = cqdg_poly + z(:,1)*g(1)

  end select
end function cqdg_poly


!    Name : cqsq_poly
! Purpose : This function computes a matrix polynomial that consists
!         : only of even order terms.
!   Input : - 2*n is the order of the polynomial.
!         : - g(i) stores the coefficient of the 2*i-th power term.
!         : - z(*,i) is the (2*i)-th power of a matrix.
!  Output : The value of the polynomial.
!    Note : This routine treats a complex, quadruple precision, square matrix.

pure function cqsq_poly(n,g,z)
  integer,          intent(in) :: n
  real   (kind=qp), intent(in) :: g(0:)
  complex(kind=qp), intent(in) :: z(:,:,0:)
  complex(kind=qp)             :: cqsq_poly(size(z,1),size(z,2))
  complex(kind=qp)             :: temporary(size(z,1),size(z,2))

  cqsq_poly = cmplx(0.0_qp, 0.0_qp, qp)

  ! an even polynomial of 2n-th order

  select case (n)                ! z(*,i) is (2*i)-th power of the matrix

    case (8)                     ! a 16-th order polynomial

      cqsq_poly =             z(:,:,0)*g(0)
      cqsq_poly = cqsq_poly + z(:,:,1)*g(1)
      cqsq_poly = cqsq_poly + z(:,:,2)*g(2)
      cqsq_poly = cqsq_poly + z(:,:,3)*g(3)

      temporary =             z(:,:,0)*g(4)
      temporary = temporary + z(:,:,1)*g(5)
      temporary = temporary + z(:,:,2)*g(6)
      temporary = temporary + z(:,:,3)*g(7)
      temporary = temporary + z(:,:,4)*g(8)
      temporary = temporary .sqtimes. z(:,:,4)

      cqsq_poly = cqsq_poly + temporary

    case (6)                     ! a 12-th order polynomial

      cqsq_poly =             z(:,:,0)*g(0)
      cqsq_poly = cqsq_poly + z(:,:,1)*g(1)
      cqsq_poly = cqsq_poly + z(:,:,2)*g(2)

      temporary =             z(:,:,0)*g(3)
      temporary = temporary + z(:,:,1)*g(4)
      temporary = temporary + z(:,:,2)*g(5)
      temporary = temporary + z(:,:,3)*g(6)
      temporary = temporary .sqtimes. z(:,:,3)

      cqsq_poly = cqsq_poly + temporary

    case (5)                     ! a 10-th order polynomial

      cqsq_poly =             z(:,:,0)*g(0)
      cqsq_poly = cqsq_poly + z(:,:,1)*g(1)
      cqsq_poly = cqsq_poly + z(:,:,2)*g(2)

      temporary =             z(:,:,0)*g(3)
      temporary = temporary + z(:,:,1)*g(4)
      temporary = temporary + z(:,:,2)*g(5)
      temporary = temporary .sqtimes. z(:,:,3)

      cqsq_poly = cqsq_poly + temporary

    case (4)                     ! an 8-th order polynomial

      cqsq_poly =             z(:,:,0)*g(0)
      cqsq_poly = cqsq_poly + z(:,:,1)*g(1)

      temporary =             z(:,:,0)*g(2)
      temporary = temporary + z(:,:,1)*g(3)
      temporary = temporary + z(:,:,2)*g(4)
      temporary = temporary .sqtimes. z(:,:,2)

      cqsq_poly = cqsq_poly + temporary

    case (3)                     ! a 6-th order polynomial

      cqsq_poly =             z(:,:,0)*g(0)
      cqsq_poly = cqsq_poly + z(:,:,1)*g(1)
      cqsq_poly = cqsq_poly + z(:,:,2)*g(2)
      cqsq_poly = cqsq_poly + z(:,:,3)*g(3)

    case (2)                     ! a 4-th order polynomial

      cqsq_poly =             z(:,:,0)*g(0)
      cqsq_poly = cqsq_poly + z(:,:,1)*g(1)
      cqsq_poly = cqsq_poly + z(:,:,2)*g(2)

    case (1)                     ! a 2-nd order polynomial

      cqsq_poly =             z(:,:,0)*g(0)
      cqsq_poly = cqsq_poly + z(:,:,1)*g(1)

  end select
end function cqsq_poly


!    Name : cqtr_poly
! Purpose : This function computes a matrix polynomial that consists
!         : only of even order terms.
!   Input : - 2*n is the order of the polynomial.
!         : - g(i) stores the coefficient of the 2*i-th power term.
!         : - z(*,i) is the (2*i)-th power of a matrix.
!  Output : The value of the polynomial.
!    Note : This routine treats a complex, quadruple precision, triangular matrix.

pure function cqtr_poly(n,g,z)
  integer,          intent(in) :: n
  real   (kind=qp), intent(in) :: g(0:)
  complex(kind=qp), intent(in) :: z(:,:,0:)
  complex(kind=qp)             :: cqtr_poly(size(z,1),size(z,2))
  complex(kind=qp)             :: temporary(size(z,1),size(z,2))

  integer :: j, pmax

  pmax = 0
  cqtr_poly = cmplx(0.0_qp, 0.0_qp, qp)

  ! an even polynomial of 2n-th order

  select case (n)                ! z(*,i) is (2*i)-th power of the matrix

    case (8)                     ! a 16-th order polynomial

      cqtr_poly =             z(:,:,0)*g(0)
      cqtr_poly = cqtr_poly + z(:,:,1)*g(1)
      cqtr_poly = cqtr_poly + z(:,:,2)*g(2)
      cqtr_poly = cqtr_poly + z(:,:,3)*g(3)

      temporary =             z(:,:,0)*g(4)
      temporary = temporary + z(:,:,1)*g(5)
      temporary = temporary + z(:,:,2)*g(6)
      temporary = temporary + z(:,:,3)*g(7)
      temporary = temporary + z(:,:,4)*g(8)
      temporary = temporary .trtimes. z(:,:,4)

      cqtr_poly = cqtr_poly + temporary
      pmax = 4

    case (6)                     ! a 12-th order polynomial

      cqtr_poly =             z(:,:,0)*g(0)
      cqtr_poly = cqtr_poly + z(:,:,1)*g(1)
      cqtr_poly = cqtr_poly + z(:,:,2)*g(2)

      temporary =             z(:,:,0)*g(3)
      temporary = temporary + z(:,:,1)*g(4)
      temporary = temporary + z(:,:,2)*g(5)
      temporary = temporary + z(:,:,3)*g(6)
      temporary = temporary .trtimes. z(:,:,3)

      cqtr_poly = cqtr_poly + temporary
      pmax = 3

    case (5)                     ! a 10-th order polynomial

      cqtr_poly =             z(:,:,0)*g(0)
      cqtr_poly = cqtr_poly + z(:,:,1)*g(1)
      cqtr_poly = cqtr_poly + z(:,:,2)*g(2)

      temporary =             z(:,:,0)*g(3)
      temporary = temporary + z(:,:,1)*g(4)
      temporary = temporary + z(:,:,2)*g(5)
      temporary = temporary .trtimes. z(:,:,3)

      cqtr_poly = cqtr_poly + temporary
      pmax = 3

    case (4)                     ! an 8-th order polynomial

      cqtr_poly =             z(:,:,0)*g(0)
      cqtr_poly = cqtr_poly + z(:,:,1)*g(1)

      temporary =             z(:,:,0)*g(2)
      temporary = temporary + z(:,:,1)*g(3)
      temporary = temporary + z(:,:,2)*g(4)
      temporary = temporary .trtimes. z(:,:,2)

      cqtr_poly = cqtr_poly + temporary
      pmax = 2

    case (3)                     ! a 6-th order polynomial

      cqtr_poly =             z(:,:,0)*g(0)
      cqtr_poly = cqtr_poly + z(:,:,1)*g(1)
      cqtr_poly = cqtr_poly + z(:,:,2)*g(2)
      cqtr_poly = cqtr_poly + z(:,:,3)*g(3)
      pmax = 3

    case (2)                     ! a 4-th order polynomial

      cqtr_poly =             z(:,:,0)*g(0)
      cqtr_poly = cqtr_poly + z(:,:,1)*g(1)
      cqtr_poly = cqtr_poly + z(:,:,2)*g(2)
      pmax = 2

    case (1)                     ! a 2-nd order polynomial

      cqtr_poly =             z(:,:,0)*g(0)
      cqtr_poly = cqtr_poly + z(:,:,1)*g(1)
      pmax = 1

  end select

  ! for the block structure to be invariant
  do j=1,size(z,1)-1
    if ( sum(abs(z(j+1,j,1:pmax))) < tiny(1.0_qp) ) then
      cqtr_poly(j+1,j)= cmplx(0.0_qp, 0.0_qp, qp)
    end if
  end do
end function cqtr_poly

#endif

end module polynomial12th

