!     Last change:  RH   22 JAN 2002   11:10 am
!     Last change:  RH   25 Apr 2001    1:37 pm
!     Last change:  RH    7 Dec 2000    1:33 pm
! Testing elementary use of Pthreads.  Print summary of tests.
! Any failures print 'failed' and cause an abort.
! This code is part of the package "A Fortran Interface to Posix Threads,"  to be
! published in ACM-TOMS.  Authors: R. Hanson, C. Breshears, and H. Gabb.
! This is test3.f90.

     MODULE global_test3A
        USE fpthrd
        IMPLICIT NONE
        INTEGER test_number, thread_counter
        LOGICAL chunk_next, ready, finished
! This is the number of working threads launched at one time.
        INTEGER, PARAMETER :: NTHREADS = 7
! This is the problem size, for a matrix-vector product example.
! The entries in the product are computed using an unrolled loop model.
! Each thread computes some of the entries of the product.
        INTEGER, PARAMETER :: MATRIX_SIZE = 127
! This is the chunk size, wherein each thread does this much of the loop.
        INTEGER, PARAMETER :: CHUNK_SIZE = 8
        TYPE(fpthrd_t) loop_thread, test_id(NTHREADS)
        TYPE(fpthrd_mutex_t) signal_mutex
        TYPE(fpthrd_cond_t) signal_unroller, signal_worker
        TYPE(fpthrd_mutexattr_t) sma, smb
        TYPE(fpthrd_attr_t) attribute

! These are the loop limits and arguments of a matrix-vector loop.
! It is part of a work crew.
        TYPE FUNCTION_ARGUMENTS
           INTEGER loop_dope(4)
           INTEGER N
           REAL, DIMENSION(:,:), POINTER :: SMATRIX
           REAL, DIMENSION(:), POINTER :: SX, SY
        END TYPE

        TYPE(function_arguments) inb(NTHREADS), inc
    END MODULE

    Module global_test3
      USE global_test3A
        INTERFACE
          SUBROUTINE TEST7(ina)
            USE global_test3A, only : function_arguments
            IMPLICIT NONE
            TYPE(function_arguments) ina
          END SUBROUTINE
          SUBROUTINE CHUNK7(INA)
            USE global_test3A, only : function_arguments
            IMPLICIT NONE
            TYPE(function_arguments), target :: INA
          END SUBROUTINE
        END INTERFACE
!           call fpthrd_create(test_id(i), NULL, chunk7, ina, status)
!        call fpthrd_create (loop_thread, NULL, test7, INC, status)

     END MODULE

     SUBROUTINE sdsdot (N, SM, SX, SY, I)
        IMPLICIT NONE
        REAL, DIMENSION(:,:), POINTER :: SM
        REAL, DIMENSION(:), POINTER :: SX, SY
        INTEGER, INTENT(IN) :: N, I

        INTEGER J
        DOUBLE PRECISION t

! Compute the dot product of two vectors.  Accumulate the results in double
! precision.  Assign the final result in single precision.
        t=0D0
        DO J=1,N
           t=t+SM(I,J)*SX(J);
        END DO
        sy(I)=t
     END SUBROUTINE

     SUBROUTINE CHUNK7(ARG_IN)
        USE global_test3, dummy => chunk7
        IMPLICIT NONE
        INTEGER :: i, local_loops, local_loope, status, summary=0
        TYPE(function_arguments), target :: arg_in
        TYPE(function_arguments), pointer :: ina

        Interface
           SUBROUTINE sdsdot (N, SM, SX, SY, I)
              IMPLICIT NONE
              REAL, DIMENSION(:,:), POINTER :: SM
              REAL, DIMENSION(:), POINTER :: SX, SY
              INTEGER, INTENT(IN) :: N, I
           END SUBROUTINE
        end interface

! Get the structure holding the arguments for this function.
        ina=>arg_in ! The variable on the left "becomes" the one on the right.
 LOOP1: DO
           call fpthrd_mutex_lock(signal_mutex, status)
           summary=summary+status

           DO while(.NOT. READY)
              call fpthrd_cond_wait(signal_worker, signal_mutex, status)
              summary=summary+status
              if(summary > 0) EXIT LOOP1
           END DO

           if(finished) THEN
              call fpthrd_mutex_unlock(signal_mutex, status)
              summary=summary+status
              EXIT LOOP1
           END IF

           local_loops=ina%loop_dope(1);local_loope=ina%loop_dope(2)
           READY=.FALSE.
           chunk_next=.TRUE.
           call fpthrd_cond_signal (signal_unroller, status)
           summary=summary+status
           call fpthrd_mutex_unlock(signal_mutex, status)
           summary=summary+status
             call ferr_abort(test_number, summary, " error with signal, lock or unlock")

! --------------------------------------------------------------------
! This is the local chunk of the loop performed by an individual worker.
           DO I=LOCAL_LOOPS, LOCAL_LOOPE
              call sdsdot(ina % n, ina%smatrix, ina % sx, ina%SY, I)
           END DO
! --------------------------------------------------------------------
        END DO LOOP1
          call ferr_abort(test_number, summary, " loop ending error with signal, lock or unlock")
     END SUBROUTINE

     SUBROUTINE TEST7(ina)
        USE global_test3, dummy => test7
        IMPLICIT NONE

        INTEGER :: nthread, chunksize, i, loops, loope, p, q, process, status, summary=0
        TYPE(function_arguments) ina
        TYPE(C_NULL) result

 ! Start up the working threads.
        call fpthrd_cond_init (signal_unroller, NULL, status)
        summary=summary+status
        call fpthrd_cond_init (signal_worker, NULL, status)
        summary=summary+status
          call ferr_abort(test_number, summary, " error with condition init")
        call ftest_comment (test_number, ". Initialized condition variables.")
! Start exercising mutex attribute manipulation functions.
        call fpthrd_mutexattr_init(sma, status)
          call ferr_abort(test_number, STATUS, " error with mutexattr init")
        call ftest_comment (test_number, ". Initialized mutex attribute.")

        call fpthrd_mutexattr_setpshared(sma, FPTHRD_PROCESS_PRIVATE, status)
          call ferr_abort(test_number, status, " error with mutexattr setpshared")

        call fpthrd_mutexattr_getpshared(sma, process, status)
          call ferr_abort(test_number, status, " error with mutexattr getpshared")

        call ftest_comment (test_number, ". Got process sharing thread.")
        if(process == FPTHRD_PROCESS_PRIVATE) THEN
           call fpthrd_mutexattr_setpshared(sma, FPTHRD_PROCESS_PRIVATE, status)
           call ferr_abort(test_number, status, " error with mutexattr setpshared")
        END IF
        call ftest_comment (test_number, ". End of mutex attribute exercises.")
! End of mutex attribute exercises.

!  May be needed in order to allow FPTHRD_mutex_setprioceiling() function to
!  operate correctly, that is, allow thread the privilege of setting
!  the priority ceiling, avoiding an EPERM error.
        smb=sma
        call fpthrd_mutexattr_setprotocol(sma, FPTHRD_PRIO_PROTECT, status)
           IF(status == ENOTSUP .or. status == ENOSYS) THEN
             sma=smb
           ELSE
             call ferr_abort(test_number, status, " error with mutexattr setprotocol")
           END IF
        call ftest_comment(test_number,". Start creating for chunk7 routine")
        call fpthrd_mutex_init(signal_mutex, sma, status)
            call ferr_abort(test_number, STATUS, " error with mutex initialization")
        ina=inc
        loops =ina % loop_dope(1);loope =ina %loop_dope(2)
        nthread=ina % loop_dope(3);chunksize=ina %loop_dope(4)
        p=loops;q=chunksize+p-1;if(q > loope)q=loope
        READY=.TRUE.
        chunk_next=.FALSE.
        finished=.FALSE.
        ina % loop_dope(1)=p;ina % loop_dope(2)=q
                DO I=1,NTHREAD
           call fpthrd_create(test_id(i), NULL, chunk7, ina, status)
             call ferr_abort(test_number, status, " creating loop worker threads")
        END DO

        test_number=test_number+1
        call ftest_comment(test_number,". Started loop worker threads and entered chunking phase")

! Chop loop into pieces of size CHUNKSIZE.  Allocate tasks to working threads.
 LOOP2: DO
           call fpthrd_mutex_lock(signal_mutex, status)
           summary=summary+status

           DO while(.NOT. chunk_next)
              call fpthrd_cond_wait(signal_unroller, signal_mutex, status)
              summary=summary+status
              if(summary > 0) EXIT LOOP2
           END DO

           if(q >= loope) EXIT LOOP2

           p=q+1;q=chunksize+p-1;if(q > loope)q=loope;READY=.TRUE.;chunk_next=.FALSE.
           ina % loop_dope(1)=p;ina % loop_dope(2)=q
           call fpthrd_cond_signal(signal_worker, status)
           summary=summary+status
           call fpthrd_mutex_unlock(signal_mutex, status)
           summary=summary+status

        END DO LOOP2

        READY=.TRUE.
        finished=.TRUE.
        call fpthrd_cond_broadcast(signal_worker, status)
        summary=summary+status
        call fpthrd_mutex_unlock(signal_mutex, status)
        summary=summary+status
          call ferr_abort (test_number, summary," lock, unlock, signal or broadcast")
        call fpthrd_mutexattr_destroy(sma, status)
        summary=summary+status
          call ferr_abort (test_number, summary," destroying mutex attribute")

        DO I=1,NTHREADS
           call fpthrd_join(test_id(i), result, status)
           summary=summary+status
        END DO
          call ferr_abort (test_number, summary," joining worker threads")
        test_number=test_number+1
        call ftest_comment(test_number,". Broadcast shutdown to worker threads and joined them.")
     END SUBROUTINE

     program fmain3
        USE global_test3
        implicit none

        REAL, POINTER :: matrix_a(:,:),vector(:),y_serial(:),y_thread(:)
        REAL errnorm, norm, temp

        integer i, j, n, status, value, priceiling
        TYPE(C_NULL) result

        Interface
           SUBROUTINE sdsdot (N, SM, SX, SY, I)
              IMPLICIT NONE
              REAL, DIMENSION(:,:), POINTER :: SM
              REAL, DIMENSION(:), POINTER :: SX, SY
              INTEGER, INTENT(IN) :: N, I
           END SUBROUTINE
        end interface

        call FPTHRD_setconcurrency(NTHREADS+1, status)

        test_number=21

! Test 7. Start a peer function that chops up a loop into threads.
! Arguments are passed by reference and are packed into a structure.
! This is a typical situation:  Use a structure to pack all arguments into one
! object.  Then pass the peer code that object.  It unpacks the structure to get
! the arguments and unrolls the loop.
        n=MATRIX_SIZE
        ALLOCATE(MATRIX_A(N,N), VECTOR(N), Y_THREAD(N), Y_SERIAL(N),STAT=status)
        call ferr_abort (test_number, status, " allocating array memory")

        inc % loop_dope(1)=1;inc%loop_dope(2)=n
        inc % loop_dope(3)=NTHREADS;inc%loop_dope(4)=CHUNK_SIZE

        inc % n=n;inc % smatrix=>matrix_a
        inc % sx=>vector;inc % sy=>y_thread
        call random_number(matrix_a)
        call random_number(vector)

        call ftest_comment (test_number, ". Launched a single startup thread.")

! Use threads to compute each row (times) vector concurrently.
        call fpthrd_create (loop_thread, NULL, test7, INC, status)
          call ferr_abort (test_number, status, "creating loop thread")

! Compute the matrix-vector product for comparison.
        DO I=1,N
           call sdsdot(n, matrix_a, vector, y_serial, I)
        END DO
! Wait for thread that unrolled and computed the loop.
        call fpthrd_join(loop_thread, result, status)
          call ferr_abort (test_number, status, "joining loop thread")

! Lock mutex to test mutex_trylock.
        call fpthrd_mutex_lock(signal_mutex, status);
          call ferr_abort (test_number, status, "locking before trymutex_lock")

! This mutex should be locked (busy) and return immediately.
        call fpthrd_mutex_trylock(signal_mutex, status)
        if(status /= EBUSY) call ferr_abort (test_number, status, "trying mutex lock when it is locked")

! Unlock mutex to clear. This may not be required.
        call fpthrd_mutex_unlock(signal_mutex, status)
          call ferr_abort (test_number, status, "unlocking after trymutex_lock")

! This mutex should be unlocked and return immediately.
                status = EBUSY
                do while (status == EBUSY)
                  call fpthrd_mutex_trylock(signal_mutex, status)
                end do
                  call ferr_abort (test_number, status, "trying mutex lock")

! Thread must have gotten lock.  It is unlocked for mutex_prioceiling tests.
                call fpthrd_mutex_unlock(signal_mutex, status)
                  call ferr_abort (test_number, status, "unlocking mutex lock after try")

! Get this mutex's priority ceiling.
        call fpthrd_mutex_getprioceiling(signal_mutex, priceiling, status)
        if(status /= ENOSYS .and. status /= ENOTSUP)&
          call ferr_abort(test_number,status," getting priority ceiling for mutex")

! Try to set the new priority to the old one.  This may gracefully fail.
        call fpthrd_mutex_setprioceiling(signal_mutex, priceiling, value, status)
        if(status /= ENOSYS .and. status /= ENOTSUP)&
          call ferr_abort(test_number,status," setting priority ceiling for mutex")

        errnorm=sum((Y_serial-y_thread)**2)
        norm=sum(y_serial**2)

! The results are correct even if they do not completely agree.
! This test will fail only with a blunder.  Relative norm errors up to about
! SQRT(EPSILON()) are allowed.
        test_number=test_number+1
        if(errnorm <= EPSILON(NORM)*norm)THEN
           call ftest_comment (test_number, ". Matrix-vector product completed; each chunk used a separate thread.")
        else
           status=1
             call ferr_abort(test_number,status," Serial and threaded matrix-vector product gave different results")
        END IF
        call skip
     END PROGRAM

