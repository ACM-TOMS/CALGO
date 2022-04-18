! This module provides random number generators of various types
! using the intrinsic subroutine, ``random_number''.
!
! This module is used only in the example programs,
! and the library itself does not require this module.
!
! The module is intended for internal use only.

module randomnumber

  use floattypes

  implicit none


  !    Name : rrandom ( a generic name )
  !   Usage : value = rrandom(s, m)
  ! Purpose : This function provides a real random number uniformly distributed
  !         : on an interval [m-s/2, m+s/2].
  !   Input : ``s'' is a real number that specifies the length of the interval.
  !         : ``m'' is a real number that specifies the center of the interval.
  !  Output : a real random number.

  interface rrandom
    module procedure rs_random                     ! in single    precision
    module procedure rw_random                     ! in double    precision
#ifdef __USE_TPREC
    module procedure rt_random                     ! in triple    precision
#endif
#ifdef __USE_QPREC
    module procedure rq_random                     ! in quadruple precision
#endif
  end interface rrandom


  !    Name : crandom ( a generic name )
  !   Usage : value = crandom(s, m)
  ! Purpose : This function provides a complex random number.
  !         : Each of the real and the imaginary part is a random number
  !         : uniformly distributed on an interval [m-s/2, m+s/2].
  !   Input : ``s'' is a real number that specifies the length of the interval.
  !         : ``m'' is a real number that specifies the center of the interval.
  !  Output : a complex random number.

  interface crandom
    module procedure cs_random                     ! in single    precision
    module procedure cw_random                     ! in double    precision
#ifdef __USE_TPREC
    module procedure ct_random                     ! in triple    precision
#endif
#ifdef __USE_QPREC
    module procedure cq_random                     ! in quadruple precision
#endif
  end interface crandom

  private
  public  :: rrandom, crandom

contains


!    Name : rs_random
!   Usage : value = rs_random(s,m)
!         : This function is invoked through the generic name, ``rrandom''.
! Purpose : This function provides a real random number in single precision
!         : uniformly distributed on an interval [m-s/2, m+s/2].
!   Input : ``s'' is a real number that specifies the length of the interval.
!         : ``m'' is a real number that specifies the center of the interval.
!  Output : a real random number in single precision.

function rs_random(s,m)
  real(kind=sp), intent(in) :: m, s
  real(kind=sp)             :: rs_random, temp(1)

  call random_number(temp(1))
  rs_random = (temp(1)-0.5_sp)*s+m
end function rs_random


!    Name : rw_random
!   Usage : value = rw_random(s,m)
!         : This function is invoked through the generic name, ``rrandom''.
! Purpose : This function provides a real random number in double precision
!         : uniformly distributed on an interval [m-s/2, m+s/2].
!   Input : ``s'' is a real number that specifies the length of the interval.
!         : ``m'' is a real number that specifies the center of the interval.
!  Output : a real random number in double precision.

function rw_random(s,m)
  real(kind=wp), intent(in) :: m, s
  real(kind=wp)             :: rw_random, temp(1)

  call random_number(temp(1))
  rw_random = (temp(1)-0.5_wp)*s+m
end function rw_random

#ifdef __USE_TPREC


!    Name : rt_random
!   Usage : value = rt_random(s,m)
!         : This function is invoked through the generic name, ``rrandom''.
! Purpose : This function provides a real random number in extended double
!         : precision uniformly distributed on an interval [m-s/2, m+s/2].
!   Input : ``s'' is a real number that specifies the length of the interval.
!         : ``m'' is a real number that specifies the center of the interval.
!  Output : a real random number in extended double precision.

function rt_random(s,m)
  real(kind=tp), intent(in) :: m, s
  real(kind=tp)             :: rt_random, temp(1)

  call random_number(temp(1))
  rt_random = (temp(1)-0.5_tp)*s+m
end function rt_random

#endif

#ifdef __USE_QPREC


!    Name : rq_random
!   Usage : value = rq_random(s,m)
!         : This function is invoked through the generic name, ``rrandom''.
! Purpose : This function provides a real random number in quadruple precision
!         : uniformly distributed on an interval [m-s/2, m+s/2].
!   Input : ``s'' is a real number that specifies the length of the interval.
!         : ``m'' is a real number that specifies the center of the interval.
!  Output : a real random number in quadruple precision.

function rq_random(s,m)
  real(kind=qp), intent(in) :: m, s
  real(kind=qp)             :: rq_random, temp(1)

  call random_number(temp(1))
  rq_random = (temp(1)-0.5_qp)*s+m
end function rq_random

#endif


!    Name : cs_random
!   Usage : value = cs_random(s,m)
!         : This function is invoked through the generic name, ``crandom''.
! Purpose : This function provides a complex random number in single precision.
!         : Each of the real and the imaginary part is a random number
!         : uniformly distributed on an interval [m-s/2, m+s/2].
!   Input : ``s'' is a real number that specifies the length of the interval.
!         : ``m'' is a real number that specifies the center of the interval.
!  Output : a complex random number in single precision.

function cs_random(s,m)
  real   (kind=sp), intent(in) :: m, s
  real   (kind=sp)             :: temp(2)
  complex(kind=sp)             :: cs_random

  call random_number(temp(1))
  call random_number(temp(2))
  cs_random = cmplx( (temp(1)-0.5_sp)*s+m, (temp(2)-0.5_sp)*s+m, sp )
end function cs_random


!    Name : cw_random
!   Usage : value = cw_random(s,m)
!         : This function is invoked through the generic name, ``crandom''.
! Purpose : This function provides a complex random number in double precision.
!         : Each of the real and the imaginary part is a random number
!         : uniformly distributed on an interval [m-s/2, m+s/2].
!   Input : ``s'' is a real number that specifies the length of the interval.
!         : ``m'' is a real number that specifies the center of the interval.
!  Output : a complex random number in double precision.

function cw_random(s,m)
  real   (kind=wp), intent(in) :: m, s
  real   (kind=wp)             :: temp(2)
  complex(kind=wp)             :: cw_random

  call random_number(temp(1))
  call random_number(temp(2))
  cw_random = cmplx( (temp(1)-0.5_wp)*s+m,(temp(2)-0.5_wp)*s+m, wp )
end function cw_random

#ifdef __USE_TPREC


!    Name : ct_random
!   Usage : value = ct_random(s,m)
!         : This function is invoked through the generic name, ``crandom''.
! Purpose : This function provides a complex random number in extended double
!         : precision. Each of the real and the imaginary part is a random
!         : number uniformly distributed on an interval [m-s/2, m+s/2].
!   Input : ``s'' is a real number that specifies the length of the interval.
!         : ``m'' is a real number that specifies the center of the interval.
!  Output : a complex random number in extended double precision.

function ct_random(s,m)
  real   (kind=tp), intent(in) :: m, s
  real   (kind=tp)             :: temp(2)
  complex(kind=tp)             :: ct_random

  call random_number(temp(1))
  call random_number(temp(2))
  ct_random = cmplx( (temp(1)-0.5_tp)*s+m,(temp(2)-0.5_tp)*s+m, tp )
end function ct_random

#endif

#ifdef __USE_QPREC


!    Name : cq_random
!   Usage : value = cq_random(s,m)
!         : This function is invoked through the generic name, ``crandom''.
! Purpose : This function provides a complex random number in quadruple
!         : precision. Each of the real and the imaginary part is a random
!         : number uniformly distributed on an interval [m-s/2, m+s/2].
!   Input : ``s'' is a real number that specifies the length of the interval.
!         : ``m'' is a real number that specifies the center of the interval.
!  Output : a complex random number in quadruple precision.

function cq_random(s,m)
  real   (kind=qp), intent(in) :: m, s
  real   (kind=qp)             :: temp(2)
  complex(kind=qp)             :: cq_random

  call random_number(temp(1))
  call random_number(temp(2))
  cq_random = cmplx( (temp(1)-0.5_qp)*s+m,(temp(2)-0.5_qp)*s+m, qp )
end function cq_random

#endif

end module randomnumber

