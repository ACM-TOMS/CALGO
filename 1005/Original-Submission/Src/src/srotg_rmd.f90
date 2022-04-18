SUBROUTINE SROTG_RMD(c, s, d, aa, ba, ca, sa, da)
  implicit none
  real    :: c, s, d, aa, ba, ca, sa, da

  ! PURPOSE
  !    Calculates the reverse mode derivative of SROTG from BLAS.
  !
  ! ARGUMENTS
  !    If SROTG was called with the arguments
  !
  !       a, b, c, s
  !
  !    then the corresponding call to SROTG_RMD should begin with the arguments
  !
  !       c, s
  !
  !    with the same values as they had on exit from SROTG. Both of these
  !    arguments will remain unchanged on exit. Note that a and b are omitted.
  !    In addition the following arguments should be provided:
  !
  !       d   (input, real scalar)
  !           the d computed by SROTG and returned in the a-parameter
  !
  !       aa  (output, real scalar)
  !           aa := the adjoint of the a supplied to SROTG
  !
  !       ba  (output, real scalar)
  !           ba := the adjoint of the b supplied to SROTG
  !
  !       ca  (input, real scalar)
  !           the adjoint of the c produced by SROTG
  !
  !       sa  (input, real scalar)
  !           the adjoint of the s produced by SROTG
  !
  !       da  (input, real scalar)
  !           the adjoint of the d returned by SROTG in the a-parameter
  !
  ! NOTES
  !    a) A sel parameter is not offered. The adjoints of a and b are always
  !       computed together; they are considered to form a pair. 
  !    b) When d = 0 the adjoints of a and b are undefined and returned as 0
  !    c) Adjoints via z which SROTG computes and returns in the b-parameter
  !       (mostly due to historical reasons) are not supported.
  !
  ! OPERATIONS
  !    BLAS: d := sigma*sqrt(a^2 + b^2)
  !          c := a/d unless d=0, then c := 1
  !          s := b/d unless d=0, then s := 0
  !          where:
  !             sigma = sign(a) if |a| > |b|
  !             sigma = sign(b) if |a| <= |b|
  !    RMD:  aa := c1 + c*d1
  !          ba := s1 + s*d1
  !          where:
  !             c1 = ca/d
  !             s1 = sa/d
  !             d1 = da - s*s1 - c*c1

  real d1, c1, s1
  
  if (d /= 0) then
    c1 = ca/d
    s1 = sa/d
    d1 = da - s*s1 - c*c1
    aa = c1 + c*d1
    ba = s1 + s*d1
  endif
end SUBROUTINE SROTG_RMD

! Note:
!   Differentiating the BLAS operation gives:
!     ba += sa/d
!     da -= sa*b/d^2 i.e. da -= sa*s/d
!     aa += ca/d
!     da -= ca*a/d^2 i.e. da -= ca*c/d
!     aa += da*sigma*a/sqrt(a^2 + b^2) i.e. aa += da*a/d i.e. aa += da*c
!     ba += da*sigma*b/sqrt(a^2 + b^2) i.e. ba += da*b/d i.e. ba += da*s

