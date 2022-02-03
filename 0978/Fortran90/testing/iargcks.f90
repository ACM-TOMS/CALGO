subroutine IARGCKS( nerrs, subnam, varnam, k, ks, nout )
!
!  Level 1 BLAS test program
!  E. Anderson
!  March 20, 2015
!
!  .. Scalar Arguments ..
   character*(*) :: subnam, varnam
   integer :: k, ks, nerrs, nout
!  ..
!
!  Purpose
!  =======
!
!  IARGCKS checks if two integer scalars are identical and reports any
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
!  K       (input) INTEGER
!          The original value of the variable VARNAM.
!
!  KS      (input) INTEGER
!          The value of the variable VARNAM after calling SUBNAM.
!
!  NOUT    (input) INTEGER
!          The unit number for output.
!
!  =====================================================================
!
   if( k /= ks ) then
      write(nout,99) subnam, varnam, k, ks
      nerrs = nerrs + 1
   end if
   return
99 format( A, ': ', A, ' was incorrectly changed from ', &
           I6, ' to ', I6 )
end subroutine
