!     Last change:  CPB  22 JAN 2002   11:10 am
!
! Testing elementary use of Fortran 90 Pthread calls.  Print summary of tests.
! Any failures print 'failed' and cause an abort.
! This code is part of the package "A Fortran Interface to Posix Threads,"  to be
! published in ACM-TOMS.  Authors R. Hanson, C. Breshears, and H. Gabb.
! This is test1.f90. Last change on 31 May 2000.

     MODULE global_test1
        USE FPTHRD
        IMPLICIT NONE
! Flag for sample output in the thread cancellation tests, TEST3().
        LOGICAL :: WANT_PRINT=.FALSE.
        INTEGER :: TEST_FREQUENCY=50000
! These are control values used in the tests.
        INTEGER :: done=0, test_number, thread_counter
! This is the number of threads launched at one time.
        INTEGER, PARAMETER :: NTHREADS=32
        LOGICAL :: STARTED=.FALSE.
        INTEGER id_no(NTHREADS), COUNTER(0:NTHREADS-1)

        TYPE(FPTHRD_t) test_id(NTHREADS)
        TYPE(FPTHRD_attr_t) :: attr
        TYPE(FPTHRD_mutex_t) :: testing_mutex = FPTHRD_MUTEX_INITIALIZER
        TYPE(FPTHRD_mutex_t) :: cancel_mutex
        TYPE(FPTHRD_mutex_t) :: once_mutex = FPTHRD_MUTEX_INITIALIZER
        TYPE(FPTHRD_once_t) :: testing_once=FPTHRD_ONCE_INIT
        TYPE(FPTHRD_cond_t) launch_done

        INTERFACE
           SUBROUTINE TEST0()
           END SUBROUTINE
           SUBROUTINE TEST1(NUMBER)
              INTEGER NUMBER
           END SUBROUTINE
           SUBROUTINE TEST2(NUMBER)
              INTEGER NUMBER
           END SUBROUTINE
           SUBROUTINE ONCE2()
           END SUBROUTINE
           SUBROUTINE TEST3(NUMBER)
              INTEGER NUMBER
           END SUBROUTINE
        END INTERFACE

     END MODULE

     SUBROUTINE TEST0()
        USE global_test1, dummy => test0
        IMPLICIT NONE
        TYPE(FPTHRD_t) localid
        INTEGER value
        call FPTHRD_self (localid)
        call FPTHRD_equal(localid, test_id(1), value)
        test_number=test_number+1
        IF(value /= 0) THEN
           CALL ftest_comment(test_number,". Startup TID is TID of this thread.")
           return
        END IF
        CALL ferr_abort (test_number, value, " matching thread IDs.")

     END SUBROUTINE TEST0

     SUBROUTINE TEST1(number)
        USE global_test1, dummy => test1
        IMPLICIT NONE
        TYPE(FPTHRD_t) localid
        INTEGER i, status, value, number

        call FPTHRD_mutex_lock(testing_mutex, status)
        call ferr_abort(test_number, status, " locking mutex before ID match loop.")
        call FPTHRD_self (localid)
        DO I=1,NUMBER
           call FPTHRD_equal(localid, test_id(I), value)
           IF(value /= 0 ) THEN
              IF(I == NTHREADS) THEN
                 call ftest_comment(test_number, ". Last stored TID is TID of last indexed thread.")
              END IF
              call FPTHRD_mutex_unlock(testing_mutex, status)
              call ferr_abort(test_number, status, " unlocking mutex during ID match loop.")
              RETURN
           END IF
        END DO
        call ferr_abort (test_number, NUMBER, "matching thread IDs.")
     END SUBROUTINE TEST1

     SUBROUTINE once2()
! This routine used to set the counter to the number of threads.
        USE global_test1, dummy => once2
        IMPLICIT NONE
        test_number=test_number+1
        call ftest_comment (test_number, ". One thread is calling 'once' function.")
        thread_counter=NTHREADS
     END SUBROUTINE once2


     SUBROUTINE TEST2 (ARG)
        USE global_test1, dummy => test2
        IMPLICIT NONE
        INTEGER ARG, VALUE, TEMPORARY
        call FPTHRD_once(testing_once, once2, value)
        IF(value /= 0) THEN
           test_number=test_number+1
           call ferr_abort (test_number, value, "launching a function, just once.")
        END IF
! An alternative to using FPTHRD_ONCE().  It is based on a mutex,
! a flag, and execution of a block of code.
 BLOCK: DO
           IF(STARTED) EXIT BLOCK ! No thread need wait after flag set.
           call FPTHRD_mutex_lock(once_mutex, value)
           call ferr_abort (test_number, value, "locking mutex for one-time code.")
           IF(STARTED) THEN
             call FPTHRD_mutex_unlock(once_mutex, value)
             call ferr_abort (test_number, value, " unlocking mutex for one-time code.")
             EXIT BLOCK ! Other threads make it through later.
           END IF
           STARTED=.TRUE. ! One thread sets the flag. Execute the one-time code.
           call ftest_comment (test_number, ". One thread is executing alternate one-time code.")
           call FPTHRD_mutex_unlock(once_mutex, value)
           call ferr_abort (test_number, value, " unlocking mutex for one-time code.")
           EXIT BLOCK
        END DO BLOCK
! Count down.  Generally must be protected by a mutex to achieve the correct value 0.
        call FPTHRD_mutex_lock(testing_mutex, value)
        call ferr_abort (test_number, value, "locking mutex for counter.")
        if(arg == NTHREADS+1) THEN
           call ftest_comment (test_number, ". Last thread index (an argument) is correctly noted.")
        END IF
        temporary=thread_counter-1
        thread_counter=temporary
! If the above and following lines are interchanged the test may fail.
        call FPTHRD_mutex_unlock(testing_mutex, value)
        call ferr_abort (test_number, value, "unlocking mutex for counter.")
     END SUBROUTINE TEST2

     SUBROUTINE TEST3 (ARG)
        USE global_test1, dummy => test3
        IMPLICIT NONE

        INTEGER arg, i, state, type, value

!   The value of the argument is the natural thread number.
!   Thread 0 will not cancel.  It completes when it is finished.
!   Threads > 0  only cancel (with testcancel) when the counter
!   has multiples of test_frequency.  Thread argument NTHREADS-1 will cancel ansynchonously.

        if(ARG == 0) THEN
           call FPTHRD_setcancelstate(FPTHRD_CANCEL_DISABLE, state, value)
           if(value /= 0)test_number=test_number+1
           call ferr_abort (test_number, value, "setting cancel state")
           if(want_print) THEN
              write(*,'(" Thread index: ",I3,2x,"ACTION: ", A)')ARG, 'FPTHRD_CANCEL_DISABLE'
           END IF

        ELSE IF(ARG < NTHREADS-1) THEN
           call FPTHRD_setcancelstate(FPTHRD_CANCEL_ENABLE, state, value)
           if(value /= 0)test_number=test_number+1
           call ferr_abort (test_number, value, "setting cancel state")
           if(want_print) THEN
              write(*,'(" Thread index: ",I3,2x,"ACTION: ", A)')ARG, 'FPTHRD_CANCEL_ENABLE'
           END IF

        END IF

        IF(ARG == NTHREADS-1) THEN
           call FPTHRD_setcanceltype (FPTHRD_CANCEL_ASYNCHRONOUS, type, value)
           if(value /= 0)test_number=test_number+1
           call ferr_abort (test_number, value, "setting cancel type")
           if(want_print) THEN
              write(*,'(" Thread index: ",I3,2x,"ACTION: ", A)')ARG, 'FPTHRD_CANCEL_ASYNCHRONOUS'
           END IF

        END IF

! Alert main thread after states and type have all been set.
        call FPTHRD_mutex_lock(cancel_mutex, value)
        call ferr_abort (test_number,value, " locking mutex")
        thread_counter=thread_counter+1
        if(thread_counter == NTHREADS) THEN
           call FPTHRD_cond_signal(launch_done, value)
           call ferr_abort(test_number,value," signalling condition")
           if(want_print) THEN
              write(*,'(" A-Thread index: ",I3,2x,"ACTION: ", A)')ARG, 'Signalling MAIN Thread'
           END IF

        END IF

        call FPTHRD_mutex_unlock(cancel_mutex, value)
        call ferr_abort (test_number, value, " unlocking mutex")

! Wait for cancellation to be sent or otherwise complete.
        DO
           IF(ARG == 0) EXIT
           counter(arg)=min(huge(1)-1,counter(arg)+1)
           if(MOD(counter(arg),test_frequency) == 0 .or. counter(arg) == huge(1)-1) THEN
              call FPTHRD_testcancel()
              if(want_print) THEN
                 write(*,'(" B3-Thread index: ",I3,2x,"ACTION: ", A)')ARG, 'Tested CANCEL and did not'
              END IF
           END IF
        END DO

! Show that disabled threads will not cancel.  They do not loop forever.
        DO
           IF(DONE > 0) THEN
              if(want_print) THEN
                 write(*,'(" B-Thread index: ",I3,2x,"ACTION: ", A)')ARG, 'SAW DONE; Exited LOOP'
              END IF

              EXIT
           END IF
           counter(arg)=min(huge(1)-1,counter(arg)+1)
           if(MOD(counter(arg),test_frequency) == 0 .or. counter(arg) == huge(1)-1) THEN
              call FPTHRD_testcancel()
              if(want_print) THEN
                 write(*,'(" B-Thread index: ",I3,2x,"ACTION: ", A)')ARG, 'Waiting for DONE; Tested CANCEL and did not'
              END IF

           END IF
        END DO


        IF(ARG == 0) THEN
           call FPTHRD_mutex_lock(cancel_mutex, value)
           call ferr_abort(test_number,value," locking mutex")
           counter(arg)=0
           thread_counter=0
           call FPTHRD_cond_signal(launch_done, value)
           call ferr_abort(test_number,value," condition signalling")

           call FPTHRD_mutex_unlock(cancel_mutex, value)
           call ferr_abort(test_number,value," unlocking mutex")

           call FPTHRD_setcancelstate(FPTHRD_CANCEL_ENABLE, state, value)
           call ferr_abort(test_number,value," re-enabling cancel state.")

! This loop spends a little time and should result in eventual cancellation.
           DO
              counter(arg)=min(huge(1)-1,counter(arg)+1)
              if(MOD(counter(arg),test_frequency) == 0 .or. counter(arg) == huge(1)-1) THEN
                 call FPTHRD_testcancel()
                 if(want_print) THEN
                    write(*,'(" E-Thread index: ",I3,2x,"ACTION: ", A)')ARG, 'Waiting for testcancel'
                 END IF
              END IF
           END DO
        END IF

     END SUBROUTINE TEST3

     Program fmain1
        USE global_test1
        IMPLICIT NONE
        INTEGER STATUS, VALUE, I, START_CONCURRENCY
        TYPE(FPTHRD_t) localid
        TYPE(C_NULL) RESULT, VALUE_RETURNED
        test_number=0

        call FPTHRD_getconcurrency(START_CONCURRENCY)
        call FPTHRD_setconcurrency(NTHREADS+1, STATUS)

! Test 0 use _create, _join, _equal, and _self to show basic activity.
! This tests launches a single thread, joins up with it, and the thread
! itself testing that its ID matches the ID obtained after launch.

        call ftest_comment (test_number, ". An FPTHRD group of _create, _join, _equal and _self, starting.")
        test_number=test_number+1
        call ftest_comment (test_number, ". Start launching a thread.")
        call FPTHRD_attr_init(attr, status)

        call ferr_abort(test_number, status," initializing thread attribute")
        call fpthrd_attr_setscope(attr,FPTHRD_SCOPE_SYSTEM, status)
        print *,"fpthrd_attr_setscope status is ",status
        if(status /= ENOSYS .and. status /= EPERM)call ferr_abort(test_number, status, " setting system scope")
        call FPTHRD_create (test_id(1), attr, test0, NULL, status)

        call ftest_comment (test_number, ". Create a thread.")
          call ferr_abort (test_number, status, "creating thread")

        call FPTHRD_join (test_id(1), NULL, status)
        test_number=test_number+1
        call ftest_comment (test_number, ". Wait joining a thread.")
        call ferr_abort (test_number, status, "joining thread")
! Test 0 is complete.
        call skip

! Test 1 use _create, _join, _equal, and _self to show concurrent activity.
! This tests launches several threads, joins up with them, and the thread
! itself tests that its ID matches the last ID obtained after launch.

        test_number=test_number+1
        call ftest_comment (test_number,". A concurrent FPTHRD group of _create, _join, _equal and _self, starting.")
        test_number=test_number+1
        call ftest_comment (test_number, ". Start launching threads.")
       
        call FPTHRD_mutex_init(testing_mutex, NULL, status)

!       testing_mutex=FPTHRD_MUTEX_INITIALIZER
        if(status /= 0) test_number=test_number+1
        call ferr_abort (test_number, status, "initializing mutex")

!   This is equivalent to:
!        testing_mutex = FPTHRD_MUTEX_INITIALIZER

!       once_mutex = FPTHRD_MUTEX_INITIALIZER
        call FPTHRD_mutex_lock(testing_mutex, status)
        call ferr_abort(test_number, status,"locking mutex")
        DO I=1,NTHREADS
           id_no(i)=i
           call FPTHRD_create (test_id(i), attr, test1, id_no(i), status)
           if(status == EAGAIN) call ftest_comment (test_number, "insufficient resources")
           call ferr_abort (test_number, status, "creating thread")
        END DO
        call FPTHRD_mutex_unlock(testing_mutex, status)

        DO I=1,NTHREADS
           call FPTHRD_join (test_id(i), NULL, status)
           call ferr_abort (test_number, status, "joining thread")
        END DO
        test_number=test_number+1
        call ftest_comment (test_number, ". Wait joining threads.")

 ! Try to join a thread already completed.  This should give a clear error.
        call FPTHRD_join(test_id(1), NULL, status)
        if(status /= ESRCH) THEN
           call ferr_abort (test_number, status, "joining thread")
        else
           test_number=test_number+1
           call ftest_comment (test_number, ". Tried to join a completed thread. Noted appropriate error.")
        end if
! Test 1 is complete.
        call skip

! Test 2 use _mutex_init, _create, _once, _join, _mutex_lock, and _mutex_unlock
!  to test well-scheduled activity. This tests launches several threads, joins
!  up with them, and the main thread tests that a counter was protected against
!  "race" conditions.

        test_number=test_number+1
        call ftest_comment (test_number,". A concurrent FPTHRD group including _once, _mutex_init, _mutex_lock, _mutex_unlock, and &
        &_mutex_destroy.")

 ! This cannot be done in a Fortran declaration, so an assignment must be used.

        DO I=1,NTHREADS
           id_no(i)=i+1
           call FPTHRD_create (test_id(i), attr, test2, id_no(i), status)
           if(status == EAGAIN) call ftest_comment (test_number, " insufficient resources")
           call ferr_abort (test_number, status, " creating thread")
        END DO

        DO I=1,NTHREADS
           call FPTHRD_join (test_id(i), NULL, status)
           if(status /= 0) test_number=test_number+1
           call ferr_abort (test_number, status, "joining thread")
        END DO

        call FPTHRD_mutex_destroy(testing_mutex,status)

        test_number=test_number+1
        if(thread_counter == 0) THEN

           call ftest_comment (test_number, ". As threads join, a mutex-protected counter was cooperatively decremented to zero.")

        else

           call ferr_abort (test_number, thread_counter, "counter was not decremented to zero")
        END IF

! Test 2 is complete.
        call skip

! Test 3: use _mutex_init, _create, _join, _mutex_lock, _mutex_unlock,
!   _cond_wait, _cond_signal, _setcancelstate, _setcanceltype, and _cancel to test
!   well-scheduled activity. This tests launches several threads, cancels
!   them, and joins up.  Some threads are temporarily left for a short time so they
!   cannot be cancelled.  After this time the thread changes its state and allows
!   cancellation.

        test_number=test_number+1
        call ftest_comment (test_number, ". A concurrent FPTHRD group including _mutex_init, _mutex_lock, _mutex_unlock, _cancel, a&
        &nd _cond_wait.")
        call ftest_comment (test_number, ". The group also includes _cond_signal,_setcancelstate, _setcanceltype, and _testcancel.")


        call FPTHRD_mutex_init(cancel_mutex, NULL, status)
        if(status /= 0) test_number=test_number+1
        call ferr_abort (test_number, status, " initializing mutex")
        call FPTHRD_cond_init(launch_done, NULL, status)
        if(status /= 0) test_number=test_number+1
        call ferr_abort (test_number, status, "initializing condition variable")

! This is equivalent to:
        launch_done = FPTHRD_COND_INITIALIZER

        thread_counter=0
        counter=0
        DO I=1,NTHREADS
           id_no(I)=I-1
           call FPTHRD_create (test_id(I), attr, test3, id_no(I), status)
           if(status == EAGAIN) call ftest_comment (test_number, "insufficient resources")
           call ferr_abort (test_number, status, " creating thread")
        END DO

        call FPTHRD_mutex_lock(cancel_mutex, status)
        call ferr_abort(test_number, status, " locking cancel mutex")


        DO WHILE(thread_counter /= NTHREADS)
           call FPTHRD_cond_wait(launch_done, cancel_mutex, status)
           call ferr_abort(test_number, status, " waiting on condition with cancel mutex")
        END DO

        call FPTHRD_mutex_unlock(cancel_mutex, status)
        call ferr_abort(test_number, status, " unlocking cancel mutex")
        test_number=test_number+1
        call ftest_comment (test_number,". Select threads are prepared for cancellation.")

! Cancel all threads.  They are waiting to be cancelled or will soon be.
        DO I=2,NTHREADS
           call FPTHRD_cancel(test_id(i), status)
           if(status /= ESRCH) call ferr_abort (test_number, status, "canceling threads")
        END DO
        if(want_print) THEN
           write(*,'(" MAIN THREAD: Cancelled all threads")')
        END IF

! This global flag results in thread 1 eventually changing to a state where it is canceled.
        DONE=1
        DO I=2,NTHREADS
           call FPTHRD_join (test_id(i), value_returned, status)
           VALUE=VALUE_RETURNED

! Expect that only the threads that could be canceled will signal with FPTHRD_CANCELED.
! Thread NTHREADS will be canceled but is not required to return with FPTHRD_CANCELED set.
           if(I < NTHREADS .and. value /= FPTHRD_CANCELED) THEN
              if(status /= 0) THEN
                 test_number=test_number+1
                 call ferr_abort (test_number, value, "joining thread but not cancelled")
              END IF
           END IF
        END DO

! Cancel first thread that was placed in a cannot cancel state, then changed to cancel.
        call FPTHRD_cancel(test_id(1), status)
        if(status /= ESRCH) call ferr_abort (test_number, status, "canceling first thread")
        call FPTHRD_mutex_lock(cancel_mutex, status)
        call ferr_abort(test_number, status, " starting second wait")

        DO WHILE (THREAD_COUNTER /= 0)
           call FPTHRD_cond_wait(launch_done, cancel_mutex, status)
           call ferr_abort(test_number, status, " at second wait")
        END DO
        call FPTHRD_mutex_unlock(cancel_mutex, status)
        call ferr_abort(test_number, status, "after second wait")

! Special treatment for select threads.
        call FPTHRD_join(test_id(1), value_returned, status)
        value=VALUE_RETURNED
        if(value /= FPTHRD_CANCELED) THEN

           if(status /= 0) test_number=test_number+1
           call ferr_abort (test_number, value, "joining thread 1 but not cancelled")

        END IF

        test_number=test_number+1
        call ftest_comment (test_number,". All elgible threads successfully cancelled.")

! Clean up the mutex and condition variables:
        call FPTHRD_mutex_destroy(cancel_mutex, status)
        call ferr_abort(test_number, status, " destroying a mutex")
        call FPTHRD_cond_destroy(launch_done, status)
        call ferr_abort(test_number, status, " destroying a conditional")
! Reset the concurrency level to the value at the start.
        call FPTHRD_setconcurrency(START_CONCURRENCY, STATUS)

! Test 4 is a call to _detach and _exit.  The main thread exits, with no
! further execution. It is an error if code executes past _exit().
        call FPTHRD_self(localid)
        call FPTHRD_detach(localid, status)
        test_number=test_number+1
        call ftest_comment (test_number,". The main thread calls _detach(itself) then _exit.")
        call skip

        call FPTHRD_exit(value)
        call ferr_abort (test_number, value, ", calling _exit, should not be here")
! Test 4 is complete.
     END program

