subroutine SARGCKS( nerrs, subnam, varnam, r, rs, nout )
   use LA_CONSTANTS32, only: wp
   use LA_XISNAN
!
!  Level 1 BLAS test program
!  E. Anderson
!  May 25, 2016
!
!  .. Scalar Arguments ..
   character*(*) :: subnam, varnam
   integer :: nerrs, nout
   real(wp) :: r, rs
!  ..
!
!  Purpose
!  =======
!
!  SARGCKS checks if two real scalars are identical and reports any
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
!  R       (input) REAL
!          The original value of the variable VARNAM.
!
!  RS      (input) REAL
!          The value of the variable VARNAM after calling SUBNAM.
!
!  NOUT    (input) INTEGER
!          The unit number for output.
!
!  =====================================================================
!
   if( r /= rs .and. ( .not.LA_ISNAN( r ) .or. .not.LA_ISNAN( rs ) ) ) then
      write(nout,99) subnam, varnam, r, rs
      nerrs = nerrs + 1
   end if
   return
99 format( A, ': ', A, ' was incorrectly changed from ', &
           E16.7, ' to ', E16.7 )
end subroutine
