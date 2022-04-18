    PROGRAM rot_test
! Test program for rewritten Givens rotation algorithm.
! This test calls the generic Fortran 90 routine.
! .. Use Statements ..
      USE import_kinds, ONLY : wp => double_precision
      USE givens_rotations, ONLY : rot
! ..
! .. Parameters ..
      REAL (wp), PARAMETER :: one = 1.0E0_wp
      REAL (wp), PARAMETER :: quarter = 0.25E0_wp
      REAL (wp), PARAMETER :: three = 3.0E0_wp
      REAL (wp), PARAMETER :: two = 2.0E0_wp
      REAL (wp), PARAMETER :: zero = 0.0E0_wp
      INTEGER, PARAMETER :: clts = 1, error = -1, rescaled = 2, sltc = 0
! ..
! .. Local Scalars ..
      REAL (wp) :: c, s
      INTEGER :: i, incx, incy, k
      CHARACTER (80) :: title
! ..
! .. Local Arrays ..
      REAL (wp) :: d(2), x(2), xv(10), yv(10)
! ..
! .. External Subroutines ..
      EXTERNAL outres
! ..
! .. Intrinsic Functions ..
      INTRINSIC sqrt
! ..
      WRITE (6,'(1x,a,/,a)') &
        ' Test Results for generic Fortran 90 version of the', &
        ' standard Givens transformation (double precision)'

! The tests are driven via a data file


10    READ (5,'(a)',end=20) title

      READ (5,*) x(1:2)
      READ (5,*) d(1:2)

! if k > 0 then there is update vector data
      READ (5,*) k
      IF (k>0) THEN
        READ (5,*) incx, incy

! incx<1 and incy<1  not sensible -- ignore
        IF (incx<1 .OR. incy<1) GO TO 10
        xv = zero
        yv = zero
        READ (5,*) (xv(i),i=1,1+(k-1)*incx,incx)
        READ (5,*) (yv(i),i=1,1+(k-1)*incy,incy)
      END IF

      WRITE (6,'(//">>> ",a)') title


! If d(1) is negative in the data file we require row
! removal. This is flagged to d_rot by a negative value
! of k.

      IF (d(1) < zero) THEN
        d(1) = -d(1)
	k = -k - 1
      END IF

! calculate the actual values of w1, z1 from the
! separate d and x values. This allows comparison
! of results between standard and modified methods.

      x = sqrt(d)*x
      xv(1:1+(k-1)*incx:incx) = sqrt(d(1))*xv(1:1+(k-1)*incx:incx)
      yv(1:1+(k-1)*incy:incy) = sqrt(d(2))*yv(1:1+(k-1)*incy:incy)
      WRITE (6,'(/" Input Data")')
      WRITE (6,'(" x-vector:", 2e16.8)') x(1:2)

      CALL rot(x(1),x(2),k,xv,incx,yv,incy,c,s)
      CALL outres(c,s,k,xv,yv,incx,incy)
      GO TO 10

20  END PROGRAM rot_test

    SUBROUTINE outres(c,s,k,xv,yv,incx,incy)
! Output routine for results
! .. Use Statements ..
      USE import_kinds, ONLY : wp => double_precision
! ..
! .. Parameters ..
      REAL (wp), PARAMETER :: zero = 0.0_wp
! ..
! .. Scalar Arguments ..
      REAL (wp) :: c, s
      INTEGER :: incx, incy, k
! ..
! .. Array Arguments ..
      REAL (wp) :: xv(*), yv(*)
! ..
! .. Local Scalars ..
      INTEGER :: i
! ..
! print up results
      WRITE (6,'(/)')

      WRITE (6,'(" c and s values:", 2e16.8)') c, s

! error condition is c=s=0
! print message and exit
      IF (c==zero .AND. s==zero) THEN
        WRITE (6,'(/" Error detected in d_rot"/)')
      ELSE


! k>0 signifies that data transformation has taken place

        IF (k>0) THEN
          WRITE (6,'(/"Output vectors x and y: incx = ",i3,", incy = ",i3)') &
            incx, incy
          WRITE (6,'("x-vector:")')
          WRITE (6,'(5x,4e16.8)') (xv(i),i=1,1+(k-1)*incx,incx)
          WRITE (6,'("y-vector:")')
          WRITE (6,'(5x,4e16.8)') (yv(i),i=1,1+(k-1)*incy,incy)
        END IF
      END IF

      WRITE (6,'(/"--------------------------------------------------")')

    END SUBROUTINE outres
