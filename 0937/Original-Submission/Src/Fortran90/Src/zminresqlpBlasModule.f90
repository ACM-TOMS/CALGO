!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     File zminresqlpBlasModule.f90
!
!     This file contains the following BLAS subroutines
!        zdotc, znrm2
!     required by subroutine MINRESQLP.
!
! Contributors:
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
! 24 Sep 2007: All parameters declared with correct intent
!              to avoid compiler warnings.
! 24 Oct 2007: Use real(8) instead of double precision or -r8.
! 24 May 2011: Use a module to package the BLAS subroutines. Use real(dp)
!              instead of real(8), where dp is a constant defined in
!              minresqlpDataModule and used in other program units.
! 12 Jul 2011: Created complex version zminresqlpBlasModule.f90
!              from real version minresqlpBlasModule.f90.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module zminresqlpBlasModule

  use  zminresqlpDataModule,    only : dp, ip, zero, one, zzero
  implicit none

  public   :: zdotc, znrm2

contains

!*****************************************************************************
function zdotc ( n, cx, incx, cy, incy )
!*****************************************************************************
!
!! CDOTC forms the dot product of two vectors, conjugating the first vector.
!
!     jack dongarra, linpack,  3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
!  Reference:
!
!    Lawson, Hanson, Kincaid, Krogh,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    Algorithm 539,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.

  implicit none

  complex(kind=dp)             :: zdotc
  complex(kind=dp)             :: ctemp
  complex(kind=dp), intent(in) :: cx(*)
  complex(kind=dp), intent(in) :: cy(*)
  integer(ip)                  :: i
  integer(ip),      intent(in) :: incx
  integer(ip),      intent(in) :: incy
  integer(ip)                  :: ix
  integer(ip)                  :: iy
  integer(ip),      intent(in) :: n

  ctemp = zzero
  zdotc = zzero

  if ( n <= 0 ) then
    return
  end if

  if ( incx == 1 .and. incy == 1 ) then

    do i = 1,n
      ctemp = ctemp + conjg ( cx(i) ) * cy(i)
    end do

  else

    if ( incx >= 0 ) then
      ix = 1
    else
      ix = ( -n + 1 ) * incx + 1
    end if

    if ( incy >= 0 ) then
      iy = 1
    else
      iy = ( -n + 1 ) * incy + 1
    end if

    do i = 1, n
      ctemp = ctemp + conjg ( cx(ix) ) * cy(iy)
      ix = ix + incx
      iy = iy + incy
    end do

  end if

  zdotc = ctemp
  return
end function zdotc


!*****************************************************************************
function znrm2 ( n, x, incx )
!*****************************************************************************
!
!! SCNRM2 returns the euclidean norm of a complex(kind=dp) vector.
!
!
!  Discussion:
!
!    SCNRM2 := sqrt ( sum ( conjg ( x(1:n) ) * x(1:n) ) )
!            = sqrt ( dot_product ( x(1:n), x(1:n) ) )
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, complex(kind=dp) X(*), the vector.
!
!    Input, integer INCX, the increment between successive entries of X.
!
!    Output, real(kind=dp) SCNRM2, the norm of the vector.

  implicit none

  integer(ip), intent(in)      :: incx
  integer(ip)                  :: ix
  integer(ip), intent(in)      :: n
  real(kind=dp)                :: norm
 !real(kind=dp), parameter     :: one = 1.0_dp
  real(kind=dp)                :: scale
  real(kind=dp)                :: znrm2
  real(kind=dp)                :: ssq
  real(kind=dp)                :: temp
  complex(kind=dp), intent(in) :: x(*)


  if ( n < 1 .or. incx < 1 ) then

    norm  = zero

  else

    scale = zero
    ssq = one

    do ix = 1, 1 + ( n - 1 ) * incx, incx
      if ( real(x(ix), dp) /= zero ) then
        temp = abs ( real(x(ix), dp) )
        if ( scale < temp ) then
          ssq = one + ssq * ( scale / temp )**2
          scale = temp
        else
          ssq = ssq + ( temp / scale )**2
        end if
      end if

      if ( aimag ( x(ix) ) /= zero ) then
        temp = abs ( aimag ( x(ix) ) )
        if ( scale < temp ) then
          ssq = one + ssq * ( scale / temp )**2
          scale = temp
        else
          ssq = ssq + ( temp / scale )**2
        end if

      end if

    end do

    norm  = scale * sqrt ( ssq )

  end if

  znrm2 = norm
  return
end function znrm2

end module zminresqlpBlasModule
