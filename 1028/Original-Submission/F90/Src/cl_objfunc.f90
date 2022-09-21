MODULE CL_OBJ_FUNC_MOD
USE OMP_LIB

! Parameters specifying file/application names.
CHARACTER(LEN=20), PARAMETER :: APPLE = "./a.out"          ! *** LINE 5
CHARACTER(LEN=20), PARAMETER :: DIR_NAME = "workdir"
CHARACTER(LEN=20), PARAMETER :: APP_INFILE = "input.dat"   ! *** LINE 7
CHARACTER(LEN=20), PARAMETER :: APP_OUTFILE = "output.dat" ! *** LINE 8
! Parameter specifying the number of tasks.
INTEGER, PARAMETER :: NUM_TASKS = 4                        ! *** LINE 10

! Locks for ensuring sequential consistency.
INTEGER(KIND=OMP_LOCK_KIND) :: CHECK_DIR_LOCK ! Master lock on directories.
LOGICAL :: DIR_LOCK(NUM_TASKS) ! Array of locks on individual directories.

CONTAINS

SUBROUTINE INIT_LOCKS
! Initialize the master lock on directory updates.
CALL OMP_INIT_LOCK(CHECK_DIR_LOCK)
! Initialize the array of locks on individual directories.
DIR_LOCK(:) = .FALSE.
RETURN
END SUBROUTINE INIT_LOCKS

SUBROUTINE DESTROY_LOCKS
! Destroy the master lock on directory updates.
CALL OMP_DESTROY_LOCK(CHECK_DIR_LOCK)
RETURN
END SUBROUTINE DESTROY_LOCKS

SUBROUTINE CL_OBJ_FUNC(C, V, IERR)
! Thread safe subroutine that calls a command line program implementing
! the objective function, while matching the interface of OBJ_FUNC expected
! by VTMOP.
USE REAL_PRECISION, ONLY : R8
IMPLICIT NONE
! Input parameters.
REAL(KIND=R8), INTENT(IN) :: C(:)
! Output parameters.
REAL(KIND=R8), INTENT(OUT) :: V(:)
INTEGER, INTENT(OUT) :: IERR
! Local variables.
CHARACTER(LEN=5) :: DIR_NUM ! Directory number in string form.
INTEGER :: I ! Loop indexing variable.
INTEGER :: TASK_NUM ! Task number.
REAL(KIND=R8) :: LAST_CHECK ! Time since last check.
! External subroutine to sleep.
EXTERNAL :: SLEEP

! Initialize the task number to 0.
TASK_NUM = 0
LAST_CHECK = OMP_GET_WTIME() - 0.5_R8
DO WHILE(TASK_NUM .EQ. 0) ! Loop until a lock is acquired.
   ! Wait for half a second before reacquiring the lock. This allows
   ! finishing tasks to easily update the directory locks.
   IF (OMP_GET_WTIME() - LAST_CHECK < 0.5_R8) CYCLE
   ! Acquire the master lock on directory updates.
   CALL OMP_SET_LOCK(CHECK_DIR_LOCK)
   DO I = 1, NUM_TASKS
      ! Check each individual working directory's lock.
      IF (.NOT. DIR_LOCK(I)) THEN
         DIR_LOCK(I) = .TRUE. ! Acquire the lock on directory I.
         TASK_NUM = I ! Set the task number.
         EXIT ! Break the inner loop as soon as a lock is acquired.
      END IF
   END DO
   ! Release the master lock on directory updates.
   CALL OMP_UNSET_LOCK(CHECK_DIR_LOCK)
   LAST_CHECK = OMP_GET_WTIME()
END DO
! Create a string version of TASK_NUM.
WRITE(DIR_NUM,'(I3)') TASK_NUM

! Write C(:) to APP_INFILE with list-directed formatting.
OPEN(100*TASK_NUM, FILE=TRIM(DIR_NAME) // TRIM(ADJUSTL(DIR_NUM)) // &
     "/" // TRIM(APP_INFILE))
WRITE(100*TASK_NUM, *) C(:)
CLOSE(100*TASK_NUM)

! Execute the command.
CALL EXECUTE_COMMAND_LINE( "cd " // TRIM(DIR_NAME) //      &
                           TRIM(ADJUSTL(DIR_NUM)) //       &
                           " && " // APPLE, EXITSTAT=IERR, &
                           WAIT=.TRUE. )
! Check the command for any errors, which are handled as missing values.
IF (IERR .NE. 0) THEN
   ! Release the lock on the individual directory TASK_NUM.
   CALL OMP_SET_LOCK(CHECK_DIR_LOCK)
   DIR_LOCK(TASK_NUM) = .FALSE.
   CALL OMP_UNSET_LOCK(CHECK_DIR_LOCK)
   RETURN
END IF

! Read V(:) from APP_OUTFILE with list-directed formatting.
OPEN(100*TASK_NUM, FILE=TRIM(DIR_NAME) // TRIM(ADJUSTL(DIR_NUM)) // &
     "/" // TRIM(APP_OUTFILE))
READ(100*TASK_NUM, *) V(:)
CLOSE(100*TASK_NUM)

! Release the lock on the individual directory TASK_NUM.
CALL OMP_SET_LOCK(CHECK_DIR_LOCK)
DIR_LOCK(TASK_NUM) = .FALSE.
CALL OMP_UNSET_LOCK(CHECK_DIR_LOCK)

RETURN
END SUBROUTINE CL_OBJ_FUNC

END MODULE CL_OBJ_FUNC_MOD
