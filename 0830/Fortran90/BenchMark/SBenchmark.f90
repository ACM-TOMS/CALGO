    FUNCTION elaptime(old_time)
!  Use this timing function in pairs:
!  TS = ELAPTIME(TICKS);  { Time something}; TS=ELAPTIME(TICKS)
!  The first call starts the clock by recording TICKS.
!  (Ignore TS the first call.)
!  The second call updates TICKS and computes TS, in Seconds.
!  This routine allows for 1-clock roll over.
!  This is a Fortran 90 standard routine.  It may return 0. if the
!  underlying processor has no clock.
! ..
! .. Function Return Value ..
      REAL :: elaptime
! ..
! .. Scalar Arguments ..
      REAL, INTENT (INOUT) :: old_time
! ..
! .. Local Scalars ..
      REAL :: newtime
! ..
! .. Intrinsic Functions ..
      INTRINSIC cpu_time
! ..
      CALL cpu_time(newtime)
      elaptime = newtime-old_time
      old_time=newtime
    END FUNCTION elaptime

    PROGRAM main_timer
! .. Use Statements ..
      USE import_kinds, ONLY : wp => single_precision
! ..
! .. Parameters ..
      REAL (wp), PARAMETER :: zero = 0.0_wp
      LOGICAL, PARAMETER :: details = .FALSE.
      INTEGER, PARAMETER :: mtests=4
      character (3), parameter :: month(12) = (/'Jan', 'Feb',&
	'Mar','Apr','May','Jun','Jul','Aug','Sep','Oct', &
        'Nov', 'Dec'/)
! ..
! .. Local Scalars ..
      INTEGER :: i, ierror, m, mon, n, ntimes, repeat_timings
      LOGICAL :: success
      CHARACTER (64) :: identification, file_name
      character (8):: date
      character (11):: time
! ..
! .. Local Arrays ..
      REAL (wp) :: maxtim(mtests), mintim(mtests)
      REAL (wp) :: mtime(mtests), sd(mtests)
      REAL (wp), ALLOCATABLE :: times(:,:)
      LOGICAL :: passed(mtests)
! ..
! .. External Functions ..
      LOGICAL, EXTERNAL :: test_pass
! ..
! .. External Subroutines ..
      EXTERNAL givens_test, output_run_details
! ..
! .. Intrinsic Functions ..
      INTRINSIC minloc, minval, real, sum, trim
! ..
! Setting details to true will generate detailed timings of 
! all methods for each run.
      WRITE (*,fmt='(/1x,A)') &
        ' Benchmark Givens Transformations on 2N by N random matrices:'
      WRITE (*,fmt='(1x,A)',advance='NO') ' ENTER the number of COLUMNS: >'
      READ (*,*,err=10,end=7) n
      WRITE (*,fmt='(1x,A)',advance='NO') &
        ' ENTER the number of times to test: >'
      READ (*,*,err=10) ntimes
      WRITE (*,fmt='(1x,A)',advance='NO') &
        ' ENTER platform identification comment: >'
      READ (*,fmt='(A64)',err=10) identification
      WRITE (*,fmt='(1x,A)',advance='NO') &
        ' ENTER how many times to repeat the whole test: >'
      READ (*,*) repeat_timings
      WRITE (*,fmt='(1x,A)',advance='NO') &
	' File name for results: >'
      READ (*,'(a)')file_name
      open(unit=4, file=file_name, status='unknown', iostat=ierror)
      IF (ierror/=0)THEN
	WRITE (*,'(" Error attempting to open results file: ",a)')&
	&            file_name
	STOP
      END IF
!
! Print header information
      call date_and_time(date,time)
      read(date(5:6),'(i2)')mon
      WRITE (4,fmt='(15x,''Givens Transformation Benchmark Program'')')
      WRITE (4,fmt='(22x,''Single Precision Version''/)')
      WRITE (4,fmt='(15x,''Run started: '',a2,'':'',a2,'' on '', &
                     &a2,'' '',a3,'' '',a4)') time(1:2), time(3:4), &
		     &date(7:8), month(mon), date(1:4)

      m = 2*n

      ALLOCATE (times(mtests,repeat_timings),STAT=ierror)
      IF (ierror>0) THEN
        WRITE (4,'(" Error attempting to allocate timing array")')
        STOP
      END IF
      success = .TRUE.
      DO i = 1, repeat_timings
        CALL givens_test(n,ntimes,times(1:mtests,i),passed)
        success = success .AND. test_pass(passed)
        IF (details) THEN
          CALL output_run_details(times(1:mtests,i),success)
        END IF
      END DO
! Compute average time
      mtime = zero
      DO i = 1, repeat_timings
        mtime = mtime + times(:,i)
      END DO
      mtime = mtime/real(repeat_timings,wp)
! and standard deviations
      sd = zero
      DO i = 1, repeat_timings
        sd = sd + (times(:,i)-mtime)**2
      END DO
      DO i = 1, mtests
        maxtim(i) = maxval(times(i,:))
        mintim(i) = minval(times(i,:))
      ENDDO

! Output synopsis of runs
      i = sum(minloc(mtime))
      WRITE (4,'(15x,"Platform id: ",a/)') trim(identification)
      WRITE (4,'(15x,"Least squares problem size: ",i4," by ",i4)') m, n
      WRITE (4,'(15x,"Problem solved ",i5," times per run")')ntimes
      WRITE (4,'(15x,"Timings averaged over ",i3," complete runs"/)')&
	     &repeat_timings
      WRITE (4,fmt='(25x,"Timing Details")')
      WRITE (4,fmt='(25x,"--------------"/)')
      WRITE (4,fmt='(28x,"New MG",4x,"New SG",3x,"Orig MG",3x,"Orig SG")')

100   format(a22,":",8f10.2)
      WRITE (4,100)"Average Times", mtime
      WRITE (4,100)"Standard Deviations", sd
      WRITE (4,100)"Normalized", mtime/minval(mtime)
      write (4,'()')
      WRITE (4,100)"Max Times", maxtim
      WRITE (4,100)"Min Times", mintim
      WRITE (4,100)"Normalized Min Times", mintim/minval(mintim)

      WRITE (4,fmt='(/5x,A)',advance='NO') &
             & " *** Optimal timings obtained using "
      SELECT CASE (i)
      CASE (1)
        WRITE (4,'("new MG")')
      CASE (2)
        WRITE (4,'("new SG")')
      CASE (3)
        WRITE (4,'("original Blas-1 MG")')
      CASE (4)
        WRITE (4,'("original Blas-1 SG")')
      END SELECT
      WRITE (4,'(/)')

      DEALLOCATE (times)
7     STOP

10    CONTINUE
      WRITE (4,'(" Error reading input")')

    END PROGRAM main_timer

    SUBROUTINE compute_sol(x,a)
! .. Use Statements ..
      USE import_kinds, ONLY : wp => single_precision
! ..
! .. Array Arguments ..
      REAL (wp), INTENT (IN) :: a(:,:)
      REAL (wp), INTENT (INOUT) :: x(:)
! ..
! .. Local Scalars ..
      INTEGER :: i, j, n
! ..
! .. Intrinsic Functions ..
      INTRINSIC size
! ..
      n = size(x)

      DO j = n, 1, -1
        x(j) = x(j)/a(j,j)
        DO i = 1, j - 1
          x(i) = x(i) - a(i,j)*x(j)
        END DO
      END DO

    END SUBROUTINE compute_sol

    FUNCTION test_pass(passed) RESULT (success)
! .. Function Return Value ..
      LOGICAL :: success
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: maxval = 4
! ..
! .. Array Arguments ..
      LOGICAL, INTENT (IN) :: passed(maxval)
! ..
! .. Local Scalars ..
      INTEGER :: i
      LOGICAL :: ltemp
! ..
      success = .TRUE.
      DO i = 1, maxval, 2
        ltemp = passed(i) .AND. passed(i+1)
        success = success .AND. ltemp
        IF ( .NOT. ltemp) THEN
          SELECT CASE ((i+1)/2)
          CASE (1)
            WRITE (4,*) 'Error in test: new MG then SG'
          CASE (2)
            WRITE (4,*) 'Error in test: original Blas-1-1 MG then SG'
          CASE DEFAULT
            WRITE (4,*) 'Illegal value of i ', i, ' in test_pass'
          END SELECT
        END IF
      END DO
    END FUNCTION test_pass

    SUBROUTINE output_run_details(times,success)
! .. Use Statements ..
      USE import_kinds, ONLY : wp => single_precision
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: maxval = 4
! ..
! .. Scalar Arguments ..
      LOGICAL :: success
! ..
! .. Array Arguments ..
      REAL (wp) :: times(maxval)
! ..
! .. Local Scalars ..
      INTEGER :: i
! ..
      IF (success) THEN
        WRITE (4,'(" All tests successfully ran")')
      ELSE
        WRITE (4,'(" IGNORE TIMINGS: One or more tests failed")')
      END IF

      DO i = 1, maxval, 2
        SELECT CASE ((i+1)/2)
        CASE (1)
          WRITE (4,'(1x,A,1x,2F7.2)') 'Timings: new MG then SG:' &
            , times(i), times(i+1)
        CASE (2)
          WRITE (4,'(1x,A,1x,2F7.2)') &
            &'Timings: original Blas-1-1 MG then SG:',&
	    &times(i), times(i+1)
        CASE DEFAULT
          WRITE (4,*) 'Illegal value of i ', i, ' in test_pass'
        END SELECT
      END DO

    END SUBROUTINE output_run_details


    SUBROUTINE givens_test(n,ntimes,times,passed)
! .. Use Statements ..
      USE import_kinds, ONLY : wp => single_precision
      USE givens_rotations, ONLY : s_rot, s_rotm
! ..
! .. Parameters ..
      REAL (wp), PARAMETER :: one = 1.0_wp
! ..
! .. Non-Generic Interface Blocks ..
      INTERFACE
        FUNCTION elaptime(old_time)
! ..
! .. Function Return Value ..
          REAL :: elaptime
! ..
! .. Scalar Arguments ..
          REAL, INTENT (INOUT) :: old_time
! ..
        END FUNCTION elaptime
      END INTERFACE
      INTERFACE
        SUBROUTINE compute_sol(x,a)
! .. Use Statements ..
          USE import_kinds, ONLY : wp => single_precision
! ..
! .. Array Arguments ..
          REAL (wp), INTENT (IN) :: a(:,:)
          REAL (wp), INTENT (INOUT) :: x(:)
! ..
        END SUBROUTINE compute_sol
      END INTERFACE
! ..
! .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: n, ntimes
      INTEGER, PARAMETER:: mtests=4
! ..
! .. Array Arguments ..
      REAL (wp), INTENT (OUT) :: times(mtests)
      LOGICAL, INTENT (OUT) :: passed(mtests)
! ..
! .. Local Scalars ..
      REAL (wp) :: c, epsval, s, tmg, tsg
      REAL :: old_time
      INTEGER :: i, ierror, itest, j, k, m
! ..
! .. Local Arrays ..
      REAL (wp), ALLOCATABLE :: a(:,:), b(:,:), d(:), x(:), y(:)
      REAL (wp) :: param(5)
! ..
! .. External Subroutines ..
      EXTERNAL srot, srotg, srotm, srotmg
! ..
! .. Intrinsic Functions ..
      INTRINSIC abs, epsilon, maxval, random_number, sqrt, transpose
! ..
      epsval = sqrt(epsilon(one))
      m = 2*n
      itest = 1
      old_time = 0.0

! Time using unit strides:
      ALLOCATE (a(n+1,m),b(n+1,m),d(m),x(n),y(n),STAT=ierror)
      IF (ierror>0) THEN
        WRITE (4,'(" Error attempting to allocate arrays for least squares")')
        STOP
      END IF

      CALL random_number(a)

! Construct and apply standard transformations:
      tsg = elaptime(old_time)
      DO k = 1, ntimes
        b = a
        DO j = 1, n
          DO i = j + 1, m
            CALL s_rot(b(j,j),b(j,i),n-j+1,b(j+1,j),1,b(j+1,i),1,c,s)
          END DO
        END DO
      END DO
      tsg = elaptime(old_time)
      times(itest) = tsg

! Check solutions:
      y(1:n) = b(n+1,1:n)
      CALL compute_sol(y(1:n),transpose(b(1:n,1:n)))

      tmg = elaptime(old_time)
      b = a
      DO k = 1, ntimes
        a = b
        d = one
        DO j = 1, n
          DO i = j + 1, m
! Construct and apply modified transformation:
            CALL s_rotm(d(j),d(i),a(j,j),a(j,i),n-j+1,a(j+1,j),1,a(j+1,i),1, &
              param)
          END DO
        END DO
      END DO
      tmg = elaptime(old_time)
      times(itest+1) = tmg
! Check solutions:
      x(1:n) = a(n+1,1:n)
      CALL compute_sol(x(1:n),transpose(a(1:n,1:n)))
      y = (x-y)/maxval(abs(y))

      passed(itest) = maxval(abs(y)) <= epsval
      passed(itest+1) = passed(itest)
      itest = itest + 2
! Obtain times with Level-1 BLAS, 2 calls per elimination.
! Construct and apply standard transformations:
      a = b
      tsg = elaptime(old_time)
      DO k = 1, ntimes
        b = a
        DO j = 1, n
          DO i = j + 1, m
            CALL srotg(b(j,j),b(j,i),c,s)
            CALL srot(n-j+1,b(j+1,j),1,b(j+1,i),1,c,s)
          END DO
        END DO
      END DO
      tsg = elaptime(old_time)
      times(itest) = tsg

! Check solutions:
      y(1:n) = b(n+1,1:n)
      CALL compute_sol(y(1:n),transpose(b(1:n,1:n)))

      tmg = elaptime(old_time)
      b = a
      DO k = 1, ntimes
        a = b
        d = one
        DO j = 1, n
          DO i = j + 1, m
! Construct and apply modified transformation:
            CALL srotmg(d(j),d(i),a(j,j),a(j,i),param)
            CALL srotm(n-j+1,a(j+1,j),1,a(j+1,i),1,param)
          END DO
        END DO
      END DO
      tmg = elaptime(old_time)
      times(itest+1) = tmg
! Check solutions:
      x(1:n) = a(n+1,1:n)
      CALL compute_sol(x(1:n),transpose(a(1:n,1:n)))
      y = (x-y)/maxval(abs(y))

      passed(itest) = maxval(abs(y)) <= epsval
      passed(itest+1) = passed(itest)
      itest = itest + 2

      DEALLOCATE (a,b,d,x,y, stat=ierror)

    END SUBROUTINE givens_test
