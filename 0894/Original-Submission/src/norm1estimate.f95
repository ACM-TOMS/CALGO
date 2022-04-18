! An interface to Hager--Higham 1-norm estimator implemented in LAPACK.
!
! Reference:
!
!   W. W. Hager,
!   Condition estimates,
!   SIAM Journal on Scientific and Statistical Computations,
!   Volume 5, pp. 311-316, 1984.
!
!   Nicholas J. Higham,
!   FORTRAN codes for estimating the one-norm of a real or complex matrix,
!   with applications to condition estimation,
!   ACM Transactions on Mathematical Software,
!   Volume. 14, Number 4, pp. 381-396, 1988.
!
! This module is intended for internal use only.

module norm1estimate

  use floattypes
  use matrixpwrtag
  use jspsylvester

  implicit none

  !    Name : sepi1estf ( a generic name )     
  !   Usage : sepinv = sepi1estf(A,B)
  ! Purpose : This function estimates the 1-norm of the inverse Sylvester
  !         : operator, (I (x) A - B^T (x) I)^{-1}, by Hager--Higham estimator.
  !   Input : A and B are upper quasi-triangular matrices.
  !  Output : The return value is the estimated 1-norm of
  !         : the inverse Sylvester operator.  

  interface sepi1estf
    module procedure rs_sepi1estf              ! Real    in Single    precision
    module procedure cs_sepi1estf              ! Complex in Single    precision
    module procedure rw_sepi1estf              ! Real    in Double    precision
    module procedure cw_sepi1estf              ! Complex in Double    precision
#ifdef __USE_TPREC
    module procedure rt_sepi1estf              ! Real    in Triple    precision
    module procedure ct_sepi1estf              ! Complex in Triple    precision
#endif
#ifdef __USE_QPREC
    module procedure rq_sepi1estf              ! Real    in Quadruple precision
    module procedure cq_sepi1estf              ! Complex in Quadruple precision
#endif
  end interface sepi1estf

  private
  public :: sepi1estf

  external :: slacon, dlacon, clacon, zlacon   ! Hager--Higham 1-norm estimator

contains


!    Name : rs_sepi1estf
! Purpose : This function estimates the 1-norm of
!         : the inverse Sylvester operator.
!   Input : A and B are upper quasi-triangular matrices.
!  Output : The return value is the estimated 1-norm of
!         : the inverse Sylvester operator.
!    Note : This function treats real matrices in single precision.
!         : SLACON in LAPACK is required.

function rs_sepi1estf(a,b)
  real(kind=sp), intent(in) :: a(1:,1:), b(1:,1:)
  real(kind=sp) :: rs_sepi1estf

  real(kind=sp) :: v(size(a,1)*size(b,1)), x(size(a,1)*size(b,1)), est
  integer       :: i(size(a,1)*size(b,1))
  integer       :: k, m, n, mn, shp(2), phs(1)

  real(kind=sp) ::  c(size(a,1),size(b,1))  ! non-transposed
  real(kind=sp) ::  y(size(a,1),size(b,1))
  real(kind=sp) :: ay(size(a,1),size(b,1))
  real(kind=sp) ::  d(size(b,1),size(a,1))  ! transposed
  real(kind=sp) ::  z(size(b,1),size(a,1))
  real(kind=sp) :: az(size(b,1),size(a,1))

  v = 0.0_sp; x = 0.0_sp; i  = 0
  c = 0.0_sp; y = 0.0_sp; ay = 0.0_sp
  d = 0.0_sp; z = 0.0_sp; az = 0.0_sp

  m  = size(a,1)
  n  = size(b,1)
  mn = m * n
  k  = 0                                    ! the parameter ``KASE''

  shp(1) = m; shp(2) = n; phs(1) = mn       ! for reshaping the arrays

  do

    call slacon(mn,v,x,i,est,k)             ! LAPACK's subroutine

    if      (k == 1) then                   ! Sylvester^{-1} * vec(X)

      c =  reshape(x,shp)                   ! unvec(x)
      call axuppersylv(a,b,c,y,ay)          ! solves AY-YB=C.
      x =  reshape(y,phs)                   ! vec(y)

    else if (k == 2) then                   ! Sylvester^{-T} * vec(X)

      c =  reshape(x,shp)                   ! unvec(x)
      d = -transpose(c)
      call axuppersylv(b,a,d,z,az)          ! solves BZ-ZA=D
      y =  transpose(z)
      x =  reshape(y,phs)                   ! vec(y)

    else                                    ! A good estimation is obtained.

      exit

    end if

  end do

  rs_sepi1estf = est
end function rs_sepi1estf


!    Name : cs_sepi1estf
! Purpose : This function estimates the 1-norm of
!         : the inverse Sylvester operator.
!   Input : A and B are upper quasi-triangular matrices.
!  Output : The return value is the estimated 1-norm of
!         : the inverse Sylvester operator.
!    Note : This function treats complex matrices in single precision.
!         : CLACON in LAPACK is required.

function cs_sepi1estf(a,b)
  complex(kind=sp), intent(in) :: a(1:,1:), b(1:,1:)
  real   (kind=sp)             :: cs_sepi1estf

  complex(kind=sp) :: v(size(a,1)*size(b,1)), x(size(a,1)*size(b,1))
  real   (kind=sp) :: est
  integer          :: k, m, n, mn, shp(2), phs(1)

  complex(kind=sp) ::  c(size(a,1),size(b,1))  ! non-transposed
  complex(kind=sp) ::  y(size(a,1),size(b,1))
  complex(kind=sp) :: ay(size(a,1),size(b,1))
  complex(kind=sp) ::  d(size(b,1),size(a,1))  ! transposed
  complex(kind=sp) ::  z(size(b,1),size(a,1))
  complex(kind=sp) :: az(size(b,1),size(a,1))
  complex(kind=sp) :: zero

  zero = cmplx(0.0_sp, 0.0_sp, sp)          ! initialization
  v = zero; x = zero
  c = zero; y = zero; ay = zero
  d = zero; z = zero; az = zero

  m  = size(a,1)
  n  = size(b,1)
  mn = m * n
  k  = 0                                    ! the parameter ``KASE''

  shp(1) = m; shp(2) = n; phs(1) = mn       ! for reshaping the arrays

  do

    call clacon(mn,v,x,est,k)               ! LAPACK's subroutine

    if      (k == 1) then                   ! Sylvester^{-1} * vec(X)

      c =  reshape(x,shp)                   ! unvec(x)
      call axuppersylv(a,b,c,y,ay)          ! solves AY-YB=C.
      x =  reshape(y,phs)                   ! vec(y)

    else if (k == 2) then                   ! Sylvester^{-H} * vec(X)

      c =  reshape(x,shp)                   ! unvec(x)
      d = -conjg(transpose(c))
      call axuppersylv(b,a,d,z,az)          ! solves BZ-ZA=D
      y =  conjg(transpose(z))
      x =  reshape(y,phs)                   ! vec(y)

    else                                    ! A good estimation is obtained.

      exit

    end if

  end do

  cs_sepi1estf = est
end function cs_sepi1estf


!    Name : rw_sepi1estf
! Purpose : This function estimates the 1-norm of
!         : the inverse Sylvester operator.
!   Input : A and B are upper quasi-triangular matrices.
!  Output : The return value is the estimated 1-norm of
!         : the inverse Sylvester operator.
!    Note : This function treats real matrices in double precision.
!         : DLACON in LAPACK is required.

function rw_sepi1estf(a,b)
  real   (kind=wp), intent(in) :: a(1:,1:), b(1:,1:)
  real   (kind=wp)             :: rw_sepi1estf

  real(kind=wp) :: v(size(a,1)*size(b,1)), x(size(a,1)*size(b,1)), est
  integer       :: i(size(a,1)*size(b,1))
  integer       :: k, m, n, mn, shp(2), phs(1)

  real(kind=wp) ::  c(size(a,1),size(b,1))  ! non-transposed
  real(kind=wp) ::  y(size(a,1),size(b,1))
  real(kind=wp) :: ay(size(a,1),size(b,1))
  real(kind=wp) ::  d(size(b,1),size(a,1))  ! transposed
  real(kind=wp) ::  z(size(b,1),size(a,1))
  real(kind=wp) :: az(size(b,1),size(a,1))

  v = 0.0_wp; x = 0.0_wp; i  = 0
  c = 0.0_wp; y = 0.0_wp; ay = 0.0_wp
  d = 0.0_wp; z = 0.0_wp; az = 0.0_wp

  m  = size(a,1)
  n  = size(b,1)
  mn = m * n
  k  = 0                                    ! the parameter ``KASE''

  shp(1) = m; shp(2) = n; phs(1) = mn       ! for reshaping the arrays

  do

    call dlacon(mn,v,x,i,est,k)             ! LAPACK's subroutine

    if      (k == 1) then                   ! Sylvester^{-1} * vec(X)

      c =  reshape(x,shp)                   ! unvec(x)
      call axuppersylv(a,b,c,y,ay)          ! solves AY-YB=C.
      x =  reshape(y,phs)                   ! vec(y)

    else if (k == 2) then                   ! Sylvester^{-T} * vec(X)

      c =  reshape(x,shp)                   ! unvec(x)
      d = -transpose(c)
      call axuppersylv(b,a,d,z,az)          ! solves BZ-ZA=D
      y =  transpose(z)
      x =  reshape(y,phs)                   ! vec(y)

    else                                    ! A good estimation is obtained.

      exit

    end if

  end do

  rw_sepi1estf = est
end function rw_sepi1estf


!    Name : cw_sepi1estf
! Purpose : This function estimates the 1-norm of
!         : the inverse Sylvester operator.
!   Input : A and B are upper quasi-triangular matrices.
!  Output : The return value is the estimated 1-norm of
!         : the inverse Sylvester operator.
!    Note : This function treats complex matrices in double precision.
!         : ZLACON in LAPACK is required.

function cw_sepi1estf(a,b)
  complex(kind=wp), intent(in) :: a(1:,1:), b(1:,1:)
  real   (kind=wp)             :: cw_sepi1estf

  complex(kind=wp) :: v(size(a,1)*size(b,1)), x(size(a,1)*size(b,1))
  real   (kind=wp) :: est
  integer          :: k, m, n, mn, shp(2), phs(1)

  complex(kind=wp) ::  c(size(a,1),size(b,1))  ! non-transposed
  complex(kind=wp) ::  y(size(a,1),size(b,1))
  complex(kind=wp) :: ay(size(a,1),size(b,1))
  complex(kind=wp) ::  d(size(b,1),size(a,1))  ! transposed
  complex(kind=wp) ::  z(size(b,1),size(a,1))
  complex(kind=wp) :: az(size(b,1),size(a,1))
  complex(kind=wp) :: zero

  zero = cmplx(0.0_wp, 0.0_wp, wp)          ! initialization
  v = zero; x = zero
  c = zero; y = zero; ay = zero
  d = zero; z = zero; az = zero

  m  = size(a,1)
  n  = size(b,1)
  mn = m * n
  k  = 0                                    ! the parameter ``KASE''

  shp(1) = m; shp(2) = n; phs(1) = mn       ! for reshaping the arrays

  do

    call zlacon(mn,v,x,est,k)               ! LAPACK's subroutine

    if      (k == 1) then                   ! Sylvester^{-1} * vec(X)

      c =  reshape(x,shp)                   ! unvec(x)
      call axuppersylv(a,b,c,y,ay)          ! solves AY-YB=C.
      x =  reshape(y,phs)                   ! vec(y)

    else if (k == 2) then                   ! Sylvester^{-H} * vec(X)

      c =  reshape(x,shp)                   ! unvec(x)
      d = -conjg(transpose(c))
      call axuppersylv(b,a,d,z,az)          ! solves BZ-ZA=D
      y =  conjg(transpose(z))
      x =  reshape(y,phs)                   ! vec(y)

    else                                    ! A good estimation is obtained.

      exit

    end if

  end do

  cw_sepi1estf = est

end function cw_sepi1estf

#ifdef __USE_TPREC


!    Name : rt_sepi1estf
! Purpose : This function estimates the 1-norm of
!         : the inverse Sylvester operator.
!   Input : A and B are upper quasi-triangular matrices.
!  Output : The return value is the estimated 1-norm of
!         : the inverse Sylvester operator.
!    Note : This function treats real matrices in triple precision.
!         : This precision is not truly supported.
!         : The input matrix is truncated to a lower precision,
!         : and an estimator for the lower precision is invoked.

function rt_sepi1estf(a,b)
  real   (kind=tp), intent(in) :: a(1:,1:), b(1:,1:)
  real   (kind=tp)             :: rt_sepi1estf
  rt_sepi1estf = rw_sepi1estf(real(a,wp), real(b,wp))
end function rt_sepi1estf


!    Name : ct_sepi1estf
! Purpose : This function estimates the 1-norm of
!         : the inverse Sylvester operator.
!   Input : A and B are upper quasi-triangular matrices.
!  Output : The return value is the estimated 1-norm of
!         : the inverse Sylvester operator.
!    Note : This function treats complex matrices in triple precision.
!         : This precision is not truly supported.
!         : The input matrix is truncated to a lower precision,
!         : and an estimator for the lower precision is invoked.

function ct_sepi1estf(a,b)
  complex(kind=tp), intent(in) :: a(1:,1:), b(1:,1:)
  real   (kind=tp)             :: ct_sepi1estf
  ct_sepi1estf = cw_sepi1estf(cmplx(a,kind=wp), cmplx(b,kind=wp))
end function ct_sepi1estf

#endif
#ifdef __USE_QPREC


!    Name : rq_sepi1estf
! Purpose : This function estimates the 1-norm of
!         : the inverse Sylvester operator.
!   Input : A and B are upper quasi-triangular matrices.
!  Output : The return value is the estimated 1-norm of
!         : the inverse Sylvester operator.
!    Note : This function treats real matrices in quadruple precision.
!         : This precision is not truly supported.
!         : The input matrix is truncated to a lower precision,
!         : and an estimator for the lower precision is invoked.

function rq_sepi1estf(a,b)
  real   (kind=qp), intent(in) :: a(1:,1:), b(1:,1:)
  real   (kind=qp)             :: rq_sepi1estf
  rq_sepi1estf = rw_sepi1estf(real(a,wp), real(b,wp))
end function rq_sepi1estf


!    Name : cq_sepi1estf
! Purpose : This function estimates the 1-norm of
!         : the inverse Sylvester operator.
!   Input : A and B are upper quasi-triangular matrices.
!  Output : The return value is the estimated 1-norm of
!         : the inverse Sylvester operator.
!    Note : This function treats complex matrices in quadruple precision.
!         : This precision is not truly supported.
!         : The input matrix is truncated to a lower precision,
!         : and an estimator for the lower precision is invoked.

function cq_sepi1estf(a,b)
  complex(kind=qp), intent(in) :: a(1:,1:), b(1:,1:)
  real   (kind=qp)             :: cq_sepi1estf
  cq_sepi1estf = cw_sepi1estf(cmplx(a,kind=wp), cmplx(b,kind=wp))
end function cq_sepi1estf

#endif
end module norm1estimate

