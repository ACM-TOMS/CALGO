    PROGRAM rotmtest
! Test program for rewritten modified Givens rotation algorithm.
! The calls are to the generic Fortran 90 version.
! .. Use Statements ..
      USE import_kinds, ONLY : wp => double_precision
      USE givens_rotations, ONLY : rotm
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
      REAL (wp) :: gamsq
      INTEGER :: i, incx, incy, k
      CHARACTER (80) :: title
! ..
! .. Local Arrays ..
      REAL (wp) :: d(2), param(5), x(2), xv(10), yv(10)
! ..
! .. External Subroutines ..
      EXTERNAL outres
! ..
! .. Intrinsic Functions ..
      INTRINSIC huge, min, tiny
! ..
      WRITE (6,'(1x,a,/,a)') &
        ' Test Results for generic Fortran 90 version of the', &
        ' modified Givens transformation (double precision)'

! main loop

10    READ (5,'(a)',end=20) title
      WRITE (6,'(//">>> ",a)') title

      READ (5,*) x(1:2)
      READ (5,*) d(1:2)


! The first tests are driven via a data file

! if k > 0 then there is update vector data

      READ (5,*) k
      IF (k>0) THEN
        READ (5,*) incx, incy
        xv = zero
        yv = zero
        READ (5,*) (xv(i),i=1,1+(k-1)*incx,incx)
        READ (5,*) (yv(i),i=1,1+(k-1)*incy,incy)
      END IF

      WRITE (6,'(/" Input Data")')
      WRITE (6,'(" x-vector:", 2e16.8)') x(1:2)
      WRITE (6,'(" d-vector:", 2e16.8)') d(1:2)


! Generate the reciprocals -- we assume that it is the d
! values that are provided. Anyway in most cases these 
! start as one so there is only a risk of confusion
! on exit.

      IF (d(1)/=zero) THEN
        d(1) = one/d(1)
      END IF
      IF (d(2)/=zero) THEN
        d(2) = one/d(2)
      END IF

      CALL rotm(d(1),d(2),x(1),x(2),k,xv,incx,yv, incy,param)
      CALL outres(param,x,d,k,xv,yv,incx,incy)
      GO TO 10

! tests for rescaling and tests for large and small data values

20    CONTINUE
      gamsq = min(huge(one),one/tiny(one))*quarter

! scaling test d(1) = d(2) < 1/g^2

      x(1) = one
      x(2) = -one
      d(1) = one/(two*gamsq)
      d(2) = d(1)
      d = one/d
      k = 2
      xv(1) = -one
      xv(2) = two
      yv(1) = one
      yv(2) = -two
      title = 'Scaling test d(1) = d(2) = 1/(2*g^2)'
      WRITE (6,'(//">>> ",a)') title
      WRITE (6,'(/" Input Data")')
      WRITE (6,'(" x-vector:", 2e16.8)') x(1:2)
      WRITE (6,'(" d-vector:", 2e16.8)') d(1:2)
      CALL rotm(d(1),d(2),x(1),x(2),k,xv,1,yv,1,param)
      CALL outres(param,x,d,k,xv,yv,incx,incy)


! d(1) = d(2) = huge/four; type = anti-diagonal

      x(1) = one
      x(2) = -one
      d(1) = huge(one)*quarter
      d(2) = huge(one)*quarter
      d = one/d
      k = -1
      title = 'd(1) = d(2) = huge/four; Type = anti-diagonal'
      WRITE (6,'(//">>> ",a)') title
      WRITE (6,'(/" Input Data")')
      WRITE (6,'(" x-vector:", 2e16.8)') x(1:2)
      WRITE (6,'(" d-vector:", 2e16.8)') d(1:2)
      CALL rotm(d(1),d(2),x(1),x(2),k,xv,1,yv,1,param)
      CALL outres(param,x,d,k,xv,yv,incx,incy)

! scaling test d(1) = d(2) = two*tiny

      x(1) = one
      x(2) = -one
      d(1) = two*tiny(one)
      d(2) = two*tiny(one)
      d = one/d
      k = -1
      title = 'Scaling test d(1) = d(2) = two*tiny'
      WRITE (6,'(//">>> ",a)') title
      WRITE (6,'(/" Input Data")')
      WRITE (6,'(" x-vector:", 2e16.8)') x(1:2)
      WRITE (6,'(" d-vector:", 2e16.8)') d(1:2)
      CALL rotm(d(1),d(2),x(1),x(2),k,xv,1,yv,1,param)
      CALL outres(param,x,d,k,xv,yv,incx,incy)

! scaling test d(1) =huge/four, d(2) = two*tiny

      x(1) = one
      x(2) = -one
      d(1) = huge(one)*quarter
      d(2) = two*tiny(one)
      d = one/d
      k = -1
      title = 'Scaling test d(1) = huge/four, d(2) = two*tiny'
      WRITE (6,'(//">>> ",a)') title
      WRITE (6,'(/" Input Data")')
      WRITE (6,'(" x-vector:", 2e16.8)') x(1:2)
      WRITE (6,'(" d-vector:", 2e16.8)') d(1:2)
      CALL rotm(d(1),d(2),x(1),x(2),k,xv,1,yv,1,param)
      CALL outres(param,x,d,k,xv,yv,incx,incy)

! scaling test d(1) = two*tiny, d(2) = huge/four

      x(1) = one
      x(2) = -one
      d(1) = two*tiny(one)
      d(2) = huge(one)*quarter
      d = one/d
      k = -1
      title = 'Scaling test d(1) = two*tiny, d(2) = huge/four'
      WRITE (6,'(//">>> ",a)') title
      WRITE (6,'(/" Input Data")')
      WRITE (6,'(" x-vector:", 2e16.8)') x(1:2)
      WRITE (6,'(" d-vector:", 2e16.8)') d(1:2)
      CALL rotm(d(1),d(2),x(1),x(2),k,xv,1,yv,1,param)
      CALL outres(param,x,d,k,xv,yv,incx,incy)

! d(1) = huge/four, d(2) = one; type = anti-diagonal

      x(1) = one
      x(2) = -one
      d(1) = huge(one)*quarter
      d(2) = one
      d = one/d
      k = -1
      title = 'd(1) = huge/four, d(2) = one; type = anti-diagonal'
      WRITE (6,'(//">>> ",a)') title
      WRITE (6,'(/" Input Data")')
      WRITE (6,'(" x-vector:", 2e16.8)') x(1:2)
      WRITE (6,'(" d-vector:", 2e16.8)') d(1:2)
      CALL rotm(d(1),d(2),x(1),x(2),k,xv,1,yv,1,param)
      CALL outres(param,x,d,k,xv,yv,incx,incy)


! scaling test d(1) = -huge/four, d(2) = two*tiny

      x(1) = one
      x(2) = -one
      d(1) = -huge(one)*quarter
      d(2) = two*tiny(one)
      d = one/d
      k = -1
      title = 'Scaling test d(1) = -huge/four, d(2) = two*tiny'
      WRITE (6,'(//">>> ",a)') title
      WRITE (6,'(/" Input Data")')
      WRITE (6,'(" x-vector:", 2e16.8)') x(1:2)
      WRITE (6,'(" d-vector:", 2e16.8)') d(1:2)
      CALL rotm(d(1),d(2),x(1),x(2),k,xv,1,yv,1,param)
      CALL outres(param,x,d,k,xv,yv,incx,incy)

! d(1) = -two*tiny, d(2) = huge/four, x(1) = 0; type = anti-diagonal

      x(1) = zero
      x(2) = -one
      d(1) = -two*tiny(one)
      d(2) = huge(one)*quarter
      d = one/d
      k = -1
      title = &
        'd(1) = -two*tiny, d(2) = huge/four, x(1) = 0; Type = anti-diagonal'
      WRITE (6,'(//">>> ",a)') title
      WRITE (6,'(/" Input Data")')
      WRITE (6,'(" x-vector:", 2e16.8)') x(1:2)
      WRITE (6,'(" d-vector:", 2e16.8)') d(1:2)
      CALL rotm(d(1),d(2),x(1),x(2),k,xv,1,yv,1,param)
      CALL outres(param,x,d,k,xv,yv,incx,incy)

! d(1) = -two*tiny, d(2) = huge/four, x(2) = 0; type = diagonal

      x(1) = one
      x(2) = zero
      d(1) = -two*tiny(one)
      d(2) = huge(one)*quarter
      d = one/d
      k = -1
      title = 'd(1) = two*tiny, d(2) = -huge/four, x(2) = 0; Type = diagonal'
      WRITE (6,'(//">>> ",a)') title
      WRITE (6,'(/" Input Data")')
      WRITE (6,'(" x-vector:", 2e16.8)') x(1:2)
      WRITE (6,'(" d-vector:", 2e16.8)') d(1:2)
      CALL rotm(d(1),d(2),x(1),x(2),k,xv,1,yv,1,param)
      CALL outres(param,x,d,k,xv,yv,incx,incy)
    END PROGRAM rotmtest

    SUBROUTINE outres(param,x,d,k,xv,yv,incx,incy)
! Output routine for the returned values
! .. Use Statements ..
      USE import_kinds, ONLY : wp => double_precision
! ..
! .. Parameters ..
      REAL (wp), PARAMETER :: one = 1.0_wp
      INTEGER, PARAMETER :: clts = 1, error = -1, rescaled = 2, sltc = 0
! ..
! .. Scalar Arguments ..
      INTEGER :: incx, incy, k
! ..
! .. Array Arguments ..
      REAL (wp) :: d(2), param(5), x(2), xv(*), yv(*)
! ..
! .. Local Scalars ..
      INTEGER :: i
! ..
! .. Intrinsic Functions ..
      INTRINSIC abs, int, sqrt
! ..
! print up results
      WRITE (6,'(/)')

      SELECT CASE (int(param(1)))

      CASE (rescaled)
        WRITE (6,'(" Type -1, Rescaling has taken place")')
        WRITE (6,'(" Row 1:   ", 2e16.8)') param(2), param(4)
        WRITE (6,'(" Row 2:   ", 2e16.8/)') param(3), param(5)
        WRITE (6,'(" x-vector:", 2e16.8)') x(1:2)
        WRITE (6,'(" d-vector:", 2e16.8)') one/d(1:2)

      CASE (sltc)
        WRITE (6,'(" Type 0, H is unit diagonal")')
        WRITE (6,'(/" Off-diagonal elements:", 2e16.8)') param(3), param(4)
        WRITE (6,'( " x-vector:             ", 2e16.8)') x(1:2)
        WRITE (6,'( " d-vector:             ", 2e16.8)') one/d(1:2)

      CASE (clts)
        WRITE (6,'(" Type 1, H is unit anti-diagonal")')
        WRITE (6,'(/" Diagonal elements: ", 2e16.8)') param(2), param(5)
        WRITE (6,'( " x-vector:          ", 2e16.8)') x(1:2)
        WRITE (6,'( " d-vector:          ", 2e16.8)') one/d(1:2)

      CASE (error)
        WRITE (6,'(" Error condition detected in Rotmg call")')
        WRITE (6,'(/"--------------------------------------------------")')
        RETURN

      END SELECT

! k>0 if data has been transformed

      IF (k>0) THEN
        WRITE (6,'(/"Output vectors x and y: incx = ",i3,", incy = ",i3)') &
          incx, incy
        WRITE (6,'("x-vector:")')
        WRITE (6,'(5x,4e16.8)') (xv(i),i=1,1+(k-1)*incx,incx)
        WRITE (6,'("y-vector:")')
        WRITE (6,'(5x,4e16.8)') (yv(i),i=1,1+(k-1)*incy,incy)
        WRITE (6,'(/"Output vectors (scaled to match _rot output)")')
        xv(1:1+(k-1)*incx:incx) = xv(1:1+(k-1)*incx:incx)/sqrt(abs(d(1)))
        yv(1:1+(k-1)*incy:incy) = yv(1:1+(k-1)*incy:incy)/sqrt(abs(d(2)))
        WRITE (6,'("x-vector:")')
        WRITE (6,'(5x,4e16.8)') (xv(i),i=1,1+(k-1)*incx,incx)
        WRITE (6,'("y-vector:")')
        WRITE (6,'(5x,4e16.8)') (yv(i),i=1,1+(k-1)*incy,incy)
      END IF
      WRITE (6,'(/"--------------------------------------------------")')
    END SUBROUTINE outres
