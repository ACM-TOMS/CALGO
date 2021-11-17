subroutine SARGCKV( nerrs, subnam, varnam, n, x, y, nout )
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
   real(wp) :: x(*), y(*)
!  ..
!
!  Purpose
!  =======
!
!  SARGCKV checks if two real vectors are identical and reports any
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
!  X       (input) REAL array, dimension (N)
!          The original value of the vector VARNAM.
!
!  Y       (input) REAL array, dimension (N)
!          The value of the vector VARNAM after calling SUBNAM.
!
!  NOUT    (input) INTEGER
!          The unit number for output.
!
!  =====================================================================
!
!  .. Local Scalars ..
   integer :: i, ierrs
!  ..
   ierrs = 0
   do i = 1, n
      if( x(i) /= y(i) .and. ( .not.LA_ISNAN( x(i) ) .or. .not.LA_ISNAN( y(i) ) ) ) then
         write(nout,99) subnam, varnam, i, y(i), x(i)
         ierrs = ierrs + 1
      end if
   end do
   if( ierrs > 0 ) nerrs = nerrs + ierrs
   return
99 format( A, ': ', A, '(', I6, ') was incorrectly changed from ', &
           E16.7, ' to ', E16.7 )
end subroutine
