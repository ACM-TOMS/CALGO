    PROGRAM dp_downdate
! .. Use Statements ..
      USE import_kinds, ONLY : wp => double_precision
! ..
! .. Parameters ..
      REAL (wp), PARAMETER :: zero = 0.0_wp
! ..
! .. Local Scalars ..
      INTEGER :: i,n, ntimes
! ..
! .. Local Arrays ..
      REAL (wp) :: ulbs(2), maxulbs(2), minulbs(2), avulbs(2)
! ..
! .. External Subroutines ..
      EXTERNAL downdate_test
! ..
! .. n is the problem size: random matrices of size 2*n x n are used
! .. ntimes defines the number of matrices that are used
      n = 40
      ntimes = 10

! Print header information
      WRITE (*,fmt='(/15x,''Downdating Test for Givens Transformation ''/)')
      WRITE (*,fmt='(12x,i5," random matrices of order ",i5, " are used")') &
     &  ntimes, n
      WRITE (*,fmt='(15x,''Errors given in multiples of ulb''/)')
! initialize
      maxulbs = zero
      minulbs = huge(zero)
      avulbs = zero
      DO i=1,ntimes
        CALL downdate_test(n,ulbs)
        maxulbs(1) = max(maxulbs(1),ulbs(1))
        minulbs(1) = min(minulbs(1),ulbs(1))
        maxulbs(2) = max(maxulbs(2),ulbs(2))
        minulbs(2) = min(minulbs(2),ulbs(2))
        avulbs = avulbs + ulbs
      END DO
      WRITE (*,'(26x, ''Average'', 9x, ''Max'',11x,''Min'')')
      WRITE (*,'(3x, ''NEW SG Downdating'',3e14.4)') &
        & avulbs(1)/real(ntimes), maxulbs(1),minulbs(1) 

      WRITE (*,'(3x, ''NEW MG Downdating'',3e14.4)') &
        & avulbs(2)/real(ntimes), maxulbs(2),minulbs(2) 

      STOP

10    CONTINUE
      WRITE (4,'(" Error reading input")')

    END PROGRAM dp_downdate

    SUBROUTINE compute_sol(x,a)
! .. Use Statements ..
      USE import_kinds, ONLY : wp => double_precision
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

    SUBROUTINE downdate_test(n,ulbs)
! .. Use Statements ..
      USE import_kinds, ONLY : wp => double_precision
      USE givens_rotations, ONLY : d_rot, d_rotm
! ..
! .. Parameters ..
      REAL (wp), PARAMETER :: one = 1.0_wp
! .. Non-Generic Interface Blocks ..
      INTERFACE
        SUBROUTINE compute_sol(x,a)
! .. Use Statements ..
          USE import_kinds, ONLY : wp => double_precision
! ..
! .. Array Arguments ..
          REAL (wp), INTENT (IN) :: a(:,:)
          REAL (wp), INTENT (INOUT) :: x(:)
! ..
        END SUBROUTINE compute_sol
      END INTERFACE
! ..
! .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: n
! ..
! .. Array Arguments ..
      REAL (wp), INTENT (OUT) :: ulbs(2)
! ..
! .. Local Scalars ..
      REAL (wp) :: c, d1, s
      INTEGER :: i, ierror, j
! ..
! .. Local Arrays ..
      REAL (wp), ALLOCATABLE :: a(:,:), b(:,:), d(:), r(:), x(:), y(:)
      REAL (wp) :: param(5)
! ..
! .. External Subroutines ..
      EXTERNAL drot, drotg, drotm, drotmg
! ..
! .. Intrinsic Functions ..
      INTRINSIC abs, maxval, random_number, transpose
! ..
! Time using unit strides:
      ALLOCATE (a(n+1,n+1),b(n+1,n+1),r(n+1),d(n+1),x(n),y(n),STAT=ierror)
      IF (ierror>0) THEN
        WRITE (4,'(" Error attempting to allocate arrays for least squares")')
        STOP
      END IF

      CALL random_number(a)

! Test add then drop accuracy.  Both Modified and Standard approaches are used.
! Only unit strides are tested.

      CALL random_number(a)
      b = a

! Reduce problem to triangular form:
      DO j = 1, n
        DO i = j + 1, n + 1
          CALL d_rot(b(j,j),b(j,i),n-j+1,b(j+1,j),1,b(j+1,i),1,c,s)
        END DO
      END DO

! Base solution:
      y(1:n) = b(n+1,1:n)
      CALL compute_sol(y(1:n),transpose(b(1:n,1:n)))

! Get a new row.  Update and then downdate.
      CALL random_number(r)
      b(1:n+1,n+1) = r
      DO j = 1, n ! ADD DATA
        CALL d_rot(b(j,j),b(j,n+1),n-j+1,b(j+1,j),1,b(j+1,n+1),1,c,s)
      END DO
      b(1:n+1,n+1) = r
      DO j = 1, n ! DROP DATA
        CALL d_rot(b(j,j),b(j,n+1),-(n-j+1)-1,b(j+1,j),1,b(j+1,n+1),1,c,s)
      END DO

      x(1:n) = b(n+1,1:n)
      CALL compute_sol(x(1:n),transpose(b(1:n,1:n)))
      x = (x-y)/maxval(abs(y))
      ulbs(1) = maxval(abs(x))/epsilon(one)

! Reduce problem to triangular form:
      d = one
      DO j = 1, n
        DO i = j + 1, n + 1
          CALL d_rotm(d(j),d(i),a(j,j),a(j,i),n-j+1,a(j+1,j),1,a(j+1,i),1, &
            param)
        END DO
      END DO

! Get a new row.  Update and then downdate.
      a(1:n+1,n+1) = r
      d(n+1) = one
      DO j = 1, n ! ADD DATA
        CALL d_rotm(d(j),d(n+1),a(j,j),a(j,n+1),n-j+1,a(j+1,j),1,a(j+1,n+1),1, &
          param)
      END DO
      a(1:n+1,n+1) = r
      d(n+1) = one
      DO j = 1, n
        d1 = -d(j) ! DROP DATA
        CALL d_rotm(d1,d(n+1),a(j,j),a(j,n+1),n-j+1,a(j+1,j),1,a(j+1,n+1),1, &
          param)
        d(j) = -d1
      END DO

      x(1:n) = a(n+1,1:n)
      CALL compute_sol(x(1:n),transpose(a(1:n,1:n)))
      x = (x-y)/maxval(abs(y))
      ulbs(2) = maxval(abs(x))/epsilon(one)

      DEALLOCATE (a,b,d,r,x,y)

    END SUBROUTINE downdate_test
