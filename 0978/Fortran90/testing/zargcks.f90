subroutine ZARGCKS( nerrs, subnam, varnam, r, rs, nout )
   use LA_CONSTANTS, only: wp
   use LA_XISNAN
!
!  Level 1 BLAS test program
!  E. Anderson
!  May 25, 2016
!
!  .. Scalar Arguments ..
   character*(*) :: subnam, varnam
   integer :: nerrs, nout
   complex(wp) :: r, rs
!  ..
!
!  Purpose
!  =======
!
!  ZARGCKS checks if two complex scalars are identical and reports any
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
!          The name of the variable whose value we are checking.
!
!  R       (input) COMPLEX
!          The original value of the variable VARNAM.
!
!  RS      (input) COMPLEX
!          The value of the variable VARNAM after calling SUBNAM.
!
!  NOUT    (input) INTEGER
!          The unit number for output.
!
!  =====================================================================
!
!  .. Local Scalars ..
   real(wp) :: a1, a2, b1, b2
!  ..
   a1 = real(r)
   a2 = real(rs)
   b1 = aimag(r)
   b2 = aimag(rs)
   if( ( a1 /= a2 .and. ( .not.LA_ISNAN( a1 ) .or. .not.LA_ISNAN( a2 ) ) ) .or. &
       ( b1 /= b2 .and. ( .not.LA_ISNAN( b1 ) .or. .not.LA_ISNAN( b2 ) ) ) ) then
      write(nout,99) subnam, varnam, real(r), aimag(r), &
         real(rs), aimag(rs)
      nerrs = nerrs + 1
   end if
   return
99 format( A, ': ', A, ' was incorrectly changed from (', &
           E16.7, ',', E16.7, ') to (', E16.7, ',', E16.7, ')' )
end subroutine
