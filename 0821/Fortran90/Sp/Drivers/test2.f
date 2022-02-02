!     Last change:  CPB  22 JAN 2002   11:10 am
!     Last change:  RH   25 JAN 2002   12:00 pm
! Testing elementary use of Pthreads.  Print summary of tests.
! Any failures print 'failed' and cause an abort.
! This code is part of the package "A Fortran Interface to Posix Threads,"  to be
! published in ACM-TOMS.  Authors: R. Hanson, C. Breshears, and H. Gabb.
! This is test2.f90.

     MODULE global_test2
        USE fpthrd
        IMPLICIT NONE
! These are control values used in the tests.
        INTEGER test_number, thread_counter
! This is the number of threads launched at one time.
        INTEGER, PARAMETER :: NTHREADS=32
        INTEGER id_no(NTHREADS)

        TYPE(FPTHRD_t) test_id(NTHREADS), first_id, detach_id
        TYPE(FPTHRD_attr_t) attribute, start_attr
        TYPE(FPTHRD_mutex_t) cancel_mutex, recursion_mutex
        TYPE(FPTHRD_once_t) testing_once
        TYPE(FPTHRD_cond_t) launch_done
        TYPE(FSCHED_PARAM) param
        TYPE(fTIMESPEC) waittime
        TYPE(fsize_t) stacksize


! Declare a derived type that will carry the problem dope.
! These are the arguments of a dot product function.

        TYPE FUNCTION_ARGUMENTS
           INTEGER N
           REAL, DIMENSION(:,:), POINTER :: SMATRIX
           REAL, DIMENSION(:), POINTER :: SX
           INTEGER INCY
           REAL, dimension(:), POINTER :: SY
           INTEGER ROW_INDEX
        END TYPE

        TYPE(FUNCTION_ARGUMENTS) inb(NTHREADS), INA

        Interface
           SUBROUTINE TEST5()
           END SUBROUTINE
           SUBROUTINE TEST6(NUMBER)
             INTEGER NUMBER
           END SUBROUTINE
           SUBROUTINE ONCE6()
           END SUBROUTINE
           RECURSIVE SUBROUTINE test7(LIMITS)
              INTEGER LIMITS(2)
           END SUBROUTINE
        end interface

     END MODULE

     SUBROUTINE FSGEMV (N, SMATRIX, SX, SY, ROW_START, ROW_END)

        IMPLICIT NONE
        INTEGER N, I, J, ROW_START, ROW_END
        REAL, pointer :: SMATRIX(:,:), SX(:)
        REAL, pointer :: SY(:)
        IF(ROW_START <= 0) RETURN
        IF(ROW_END > N) RETURN
        IF(ROW_END < ROW_START) RETURN
        SY(ROW_START:ROW_END)=0E0
        DO J=1,N
           DO I=ROW_START, ROW_END
              SY(I)=SY(I)+SMATRIX(I,J)*SX(J)
           END DO
        END DO

     END SUBROUTINE

     SUBROUTINE FSDSDOT (N, SMATRIX, SX, INCY, SY, ROW_INDEX)

        IMPLICIT NONE

        INTEGER N, INCY, J, Q, ROW, ROW_INDEX
        REAL, pointer :: SMATRIX(:,:), SX(:)
        REAL, pointer :: SY(:)
        DOUBLE PRECISION T
        ROW=ROW_INDEX
        T=0D0
        Q=1

        DO J=1,N
           T=T+SMATRIX(ROW, J)*SX(Q)
           Q=Q+INCY
        END DO
        SY(ROW)=T
     END SUBROUTINE

     subroutine once6()
        use global_test2, dummy => once6
        implicit none

! Use with FPTHRD_once() to initialize any data element.
        integer status
        recursion_mutex=FPTHRD_MUTEX_INITIALIZER
        test_number=test_number+1
        call ftest_comment(test_number, ". Once initializing a mutex that is used later.")
     end subroutine


     SUBROUTINE TEST5()
        USE global_test2, dummy => test5
        IMPLICIT NONE
        integer value,  TIMES, TIMEE, RATE
! Not even a cancellation request can break out of this spin loop.
! The thread does its work, then signals the main thread it is finished.
        CALL SYSTEM_CLOCK(TIMES, COUNT_RATE=RATE)
        DO
           CALL SYSTEM_CLOCK(TIMEE)
           IF(TIMEE <= TIMES) EXIT ! Avoid a clock roll-over.
           IF((TIMEE-TIMES)/2 >= RATE) EXIT ! Loop about two seconds with this test.
        END DO

! After this thread has done its work it unlocks and signals.
        call FPTHRD_mutex_unlock(cancel_mutex, value)
          call ferr_abort(test_number, value, " unlocking mutex.")
        call FPTHRD_cond_signal(launch_done, value)
          call ferr_abort(test_number, value, " condition signal.")
     END SUBROUTINE


     subroutine test6(arg_in)
        use global_test2, dummy => test6
        implicit none
        INTERFACE
           SUBROUTINE FSDSDOT (N, SMATRIX, SX, INCY, SY, ROW_INDEX)
              INTEGER N, INCY, ROW_INDEX
              REAL, POINTER, DIMENSION(:) :: SMATRIX(:,:), SX, SY
           END SUBROUTINE
        END INTERFACE

        integer status, value, arg_in, I
  
        call FPTHRD_once (testing_once, once6, status)
          call ferr_abort(test_number, status, "starting once function for thread data key")
        I=arg_in

! Use the structure holding the arguments for this function.
        call fsdsdot (inb(I) % n, inb(I) % SMATRIX, inb(I) % SX, inb(I) % incy, inb(I) % SY, inb(I)% row_INDEX)
     end subroutine

     recursive subroutine test7(LIMITS)
        use global_test2, my_test=>test7
! Symbol is reset to a dummy. Avoids conflict.
        implicit none

        INTERFACE
           SUBROUTINE FSGEMV (N, SMATRIX, SX, SY, ROW_START, ROW_END)
              IMPLICIT NONE
              INTEGER N, ROW_START, ROW_END
              REAL, POINTER, DIMENSION(:) :: SMATRIX(:,:), SX, SY
           END SUBROUTINE
        END INTERFACE

        INTEGER :: IDEAL, K, J, LIMITS(2), LIMITS_L(2), LIMITS_R(2), status
        TYPE(FPTHRD_T) THREAD_L, THREAD_R

        call FPTHRD_mutex_lock(recursion_mutex, status)
          call ferr_abort(test_number, status, " locking mutex in test7")
        K=LIMITS(2)-LIMITS(1)+1;J=(LIMITS(1)+LIMITS(2))/2

! The problem limits are split into two equally sized groups.
        LIMITS_L=LIMITS;LIMITS_R=LIMITS
        LIMITS_L(2)=J;LIMITS_R(1)=J+1
        call FPTHRD_mutex_unlock(recursion_mutex, status)
          call ferr_abort(test_number, status, " unlocking mutex in test7")
        IDEAL=(INA%N + 7)/8
        IF(K <= IDEAL ) THEN

! This is where the work actually gets done.  The above value of IDEAL is arbitrary.
           CALL FSGEMV (INA % N, INA % SMATRIX, INA % SX, INA % SY, LIMITS(1), LIMITS(2))
        ELSE
           call FPTHRD_create(THREAD_L, attribute, my_test, limits_l, status)
             call ferr_abort(test_number, status, " recursive create-L in test7")
           call FPTHRD_create(THREAD_R, attribute, my_test, limits_r, status)
             call ferr_abort(test_number, status, " recursive create-R in test7")
           call FPTHRD_join(THREAD_L, NULL, status)
             call ferr_abort(test_number, status, " recursive join-L in test7")
           call FPTHRD_join(THREAD_R, NULL, status)
             call ferr_abort(test_number, status, " recursive join-R in test7")
        END IF

     END SUBROUTINE

     program fmain2
        USE global_test2
        implicit none

        REAL, POINTER :: matrix_a(:,:),vector(:),y_serial(:),y_thread(:)
        REAL errnorm, norm, temp

        integer create_state, i, j, n, status, value, LIMITS(2), &
          change_sec, change_nanosec, schedule_value, new_value

        test_number=14

        call FPTHRD_setconcurrency(NTHREADS+1) ! Optional argument not used.



! Test 5: Examine thread attributes.  Reset values and qualities of the threads.
!  This tests launches a single detached thread, and waits for it.

        test_number=test_number+1
        call ftest_comment (test_number, ". An FPTHRD group examining and changing attributes.")
        call FPTHRD_attr_init(attribute, status)
          call ferr_abort(test_number, status, " initializing attribute")
        call FPTHRD_attr_getstacksize(attribute, stacksize, status)
        call ferr_abort(test_number, status, " getting default stack size")
        call FPTHRD_attr_getstacksize(attribute, stacksize, status)
        call fpthrd_get_fsize(value, stacksize)
! The value FPTHRD_STACK_MIN may be set to a flag value (-1,0, etc)
! that indicates that it is at some unspecified default.
        IF(FPTHRD_STACK_MIN <= 0) THEN
          value = value+32000
        ELSE
          value = 3*FPTHRD_STACK_MIN/2
        END IF
! Reset stack size to average of the minimum and the default.
        call fpthrd_set_fsize(value, stacksize)

        call FPTHRD_attr_setstacksize(attribute, stacksize, status)
        call ferr_abort (test_number, status, "setting stack size")
        call FPTHRD_attr_getstacksize(attribute, stacksize, status)
        call fpthrd_get_fsize(new_value, stacksize)
        if(value /= new_value)& 
        call ferr_abort(test_number, 1,&
        " stacksize is not equal to 3*FPTHRD_STACK_MIN/2 or initial value.")
 
! Initialize mutex and condition variable.
        cancel_mutex=FPTHRD_MUTEX_INITIALIZER
        launch_done =FPTHRD_COND_INITIALIZER

! Create a detached thread.
        call FPTHRD_attr_init(start_attr, status)
          call ferr_abort(test_number, status," initializing attribute")

        call FPTHRD_attr_getdetachstate(start_attr, create_state, status)
        if(status == EINVAL)call ferr_abort (test_number, status, "getting default detached state")

        if(create_state /= FPTHRD_CREATE_DETACHED) THEN
           call FPTHRD_attr_setdetachstate(start_attr, FPTHRD_CREATE_DETACHED, status)
             call ferr_abort (test_number, status, "setting detached state")
        END IF


! The thread to be created unlocks the mutex so that the signal occurs.
        call FPTHRD_mutex_lock(cancel_mutex, status)
          call ferr_abort (test_number, status, "locking mutex")

        test_number=test_number+1
        call ftest_comment (test_number, ". Launch a detached thread.")
        call FPTHRD_create (detach_id, start_attr, test5, NULL, status)


        call ferr_abort (test_number, status, "creating thread")

! Wait for the detached thread to signal that it has completed.
! A little more time is allowed than what is expected.

! Assigns the current epoch plus the value on the right side.
        change_sec=3; change_nanosec=0
        Call fpthrd_set_ftimespec(change_sec, change_nanosec, waittime)
! Equivalent C code:
!   waittime.tv_sec=time(NULL)+3;
!   waittime.tv_nsec=0;
        call ftest_comment (test_number, ". Entering timed wait.")

        call FPTHRD_cond_timedwait(launch_done, cancel_mutex, waittime, status)
        if(status /= ETIMEDOUT) call ferr_abort (test_number, status, " timed wait")
        call ftest_comment (test_number, ". Left timed wait.")

        call FPTHRD_mutex_unlock(cancel_mutex, status)
          call ferr_abort (test_number, status, " unlocking mutex")

        test_number=test_number+1
        call ftest_comment (test_number, ". Meeting a detached thread after it completes.")
        if(status == 0) call ferr_abort (test_number, status, "meeting a detached thread")
        call FPTHRD_attr_destroy(start_attr, status)

! Test 5 is complete.
        call skip

! Test 6. Start a peer function that calls a library function.  Arguments
!           are packed into a derived type.
        N=NTHREADS

! No claims of randomness are made for this sequence.  It is used to generate
! a non-repeatable sequence of matrix and vector values.
        ALLOCATE(MATRIX_A(N,N), VECTOR(N), Y_SERIAL(N), Y_THREAD(N), STAT=status)
        call ferr_abort (test_number, status, "allocating array memory")
        call random_number(matrix_a)
        call random_number(vector)
        DO J=1,N
           vector(J)=2e0*vector(J)-1e0
        END DO

        testing_once=FPTHRD_ONCE_INIT
! The value FPTHRD_STACK_MIN may be set to a flag value (-1,0, etc)
! that indicates that it is at some unspecified default.
        IF(FPTHRD_STACK_MIN <= 0) THEN
          call FPTHRD_attr_getstacksize(attribute, stacksize, status)
          call fpthrd_get_fsize(value, stacksize)
          value=value*2          
        ELSE
          value = 2*FPTHRD_STACK_MIN
        END IF
        call fpthrd_set_fsize(value, stacksize)

        call FPTHRD_attr_setstacksize(attribute, stacksize, status)
        call ferr_abort (test_number, status, " setting stack size")
     
! Use threads to compute each row (times) vector concurrently.
        DO I=1,N
! This is a typical situation:  Use a structure to pack all arguments into one
! object.  Then pass the peer code that object.  It unpacks the structure to get
! the arguments.

           inb(I) % n=n;inb(I) % SMATRIX=> matrix_a;
           inb(I) % SX=>vector;inb(I) % incy=1;inb(I) % SY=>y_thread;
           inb(I) % row_index=I
           ID_NO(I)=I

! Create the set of threads for the product.
           call FPTHRD_create (test_id(i), attribute, test6, ID_NO(I), status)

             call ferr_abort (test_number, status, " creating threads")
        END DO

! Compute the matrix-vector product for comparison.
! This computation may overlap the thread computation.
        y_serial=matmul(matrix_a, vector)

        DO I=1,NTHREADS
           call FPTHRD_join(test_id(i), NULL, status)
             call ferr_abort(test_number, status, " joining a single thread")
        END DO
        errnorm=sum((y_serial-y_thread)**2)
        norm=sum(y_serial**2)

! The results are correct even if they do not completely agree.
! This test will fail only with a blunder.  Small relative errors will be allowed.
        test_number=test_number+1

        if(errnorm <= EPSILON(NORM)*norm) THEN
           call ftest_comment (test_number, ". Matrix-vector product with each entry using a separate thread.")
        else
           status=1
             call ferr_abort(test_number,status," Serial and threaded matrix-vector product gave different results")
        END IF
        call skip
! Test 6 is complete. */
        y_thread=0e0
! Use threads to compute each row (times) vector concurrently.
! This exericse uses one thread to call a routine.  The routine uses divide and
! conquer  (recursion) to reduce the problem size to one with good properties:
! (The matrix dimension K by N, where K <= N/8).

        ina % N=n;ina % SMATRIX=> matrix_a;ina % SX=>vector;ina % SY=>y_thread
        call FPTHRD_attr_setinheritsched(attribute, FPTHRD_EXPLICIT_SCHED, status)
          call ferr_abort(test_number, status, " setting inherit schedule")

        call FPTHRD_attr_setschedpolicy(attribute, FSCHED_RR, status)
! This may not be supported:
        IF(status /= ENOTSUP) call ferr_abort(test_number, status, " setting schedule policy")
        call FPTHRD_attr_setinheritsched(attribute, FPTHRD_INHERIT_SCHED, status)
! This may not be supported:
        IF(status /= ENOTSUP) call ferr_abort(test_number, status, " setting inherit schedule")
        schedule_value=1
        call fpthrd_set_fsched_param (schedule_value, param)

        call fpthrd_attr_setschedparam(attribute, param, status)
          call ferr_abort(test_number, status, " setting schedule parameter")

! Create the thread for the product.
        LIMITS=(/1,N/)
        call FPTHRD_create (first_id, attribute, test7, LIMITS, status)
          call ferr_abort (test_number, status, " creating a recursive thread")

        call FPTHRD_join(first_id, NULL, status)
          call ferr_abort(test_number, status, " joining a single recursive thread")
        call fpthrd_get_fsched_param (status, param)
! Check that schedule parameter component was communicated.
          call ferr_abort(test_number, status-schedule_value,&
            " schedule value changed and then was not recovered")
! Check results for correctness:
        errnorm=sum((y_serial-y_thread)**2)
        norm=sum(y_serial**2)
! The results are correct even if they do not completely agree.
! This test will fail only with a blunder.  Small relative errors will be allowed.
        test_number=test_number+1

        if(errnorm <= EPSILON(NORM)*norm) THEN
           call ftest_comment (test_number, ". Matrix-vector product with recursive threads.")
        else
           status=1
             call ferr_abort(test_number,status," Serial and threaded matrix-vector product gave different results")
        END IF
        call skip
! Test 7 is complete. */

     end program


