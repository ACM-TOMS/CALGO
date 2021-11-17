subroutine CARGCKV( nerrs, subnam, varnam, n, x, y, nout )
   use LA_CONSTANTS32, only: wp
   use LA_XISNAN
!
!  Level 1 BLAS test program
!  E. Anderson
!  May 25, 2016
!
!  .. Scalar Arguments ..
   character*(*) :: subnam, varnam
   integer :: n, nerrs, nout
!  ..
!  .. Array Arguments ..
   complex(wp) :: x(*), y(*)
!  ..
!
!  Purpose
!  =======
!
!  CARGCKV checks if two complex vectors are identical and reports any
!  errors.
!
!  Arguments
!  =========
!
!  NERRS   (input/output) INTEGER
!          The number of errors, incremented if any are detected in this
!          subroutine.
!
!  SUBNAM  (input) CHARACTER*(*)
!          The name of the subroutine whose arguments we are checking.
!
!  VARNAM  (input) CHARACTER*(*)
!          The name of the variable whose values we are checking.
!
!  N       (input) INTEGER
!          The length of the vector VARNAM in SUBNAM's argument list.
!
!  X       (input) COMPLEX array, dimension (N)
!          The original value of the vector VARNAM.
!
!  Y       (input) COMPLEX array, dimension (N)
!          The value of the vector VARNAM after calling SUBNAM.
!
!  NOUT    (input) INTEGER
!          The unit number for output.
!
!  =====================================================================
!
!  .. Local Scalars ..
   integer :: i, ierrs
   real(wp) :: a1, a2, b1, b2
!  ..
   ierrs = 0
   do i = 1, n
      a1 = real(x(i))
      a2 = real(y(i))
      b1 = aimag(x(i))
      b2 = aimag(y(i))
      if( ( a1 /= a2 .and. ( .not.LA_ISNAN( a1 ) .or. .not.LA_ISNAN( a2 ) ) ) .or. &
          ( b1 /= b2 .and. ( .not.LA_ISNAN( b1 ) .or. .not.LA_ISNAN( b2 ) ) ) ) then
         write(nout,99) subnam, varnam, i, real(y(i)), aimag(y(i)), &
            real(x(i)), aimag(x(i))
         ierrs = ierrs + 1
      end if
   end do
   if( ierrs > 0 ) nerrs = nerrs + ierrs
   return
99 format( A, ': ', A, '(', I6, ') was incorrectly changed from (', &
           E16.7, ',', E16.7, ') to (', E16.7, ',', E16.7, ')' )
end subroutine
